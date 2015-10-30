!---------------------------------------------------------------------------
!             AUTHOR TRISTAN SALLES -- UNIVERSITY OF SYDNEY
!---------------------------------------------------------------------------
module diffusionSolve

  use diffusionParam

  implicit none

contains

  !---------------------------------------------------------------------------
  subroutine solve_explicit_PDE

    integer::i,j,k

    real(kind=8)::sum,ax,ay,diffp,diffm,diff,th
    real(kind=8)::halfxp,halfxm,halfyp,halfym

    real(kind=8),dimension(nlith)::nlitho
    real(kind=8),dimension(nx,ny)::hrate

    ! First we solve the equation on elevation
    dt=maxstep
    do j=js(pet_id+1)+1,je(pet_id+1)-1
      do i=is(pet_id+1)+1,ie(pet_id+1)-1
        sum=0.
        do k=1,nlith
          diff=litho(k)%diffa
          if(gz(i,j)<0.0) diff=litho(k)%diffm
          diffp=litho(k)%diffa
          if(gz(i+1,j)<0.0) diffp=litho(k)%diffm
          diffm=litho(k)%diffa
          if(gz(i-1,j)<0.0) diffm=litho(k)%diffm
          halfxp=0.5*(diff*Flitho(i,j,k)+diffp*Flitho(i+1,j,k))
          halfxm=0.5*(diff*Flitho(i,j,k)+diffm*Flitho(i-1,j,k))
          diffp=litho(k)%diffa
          if(gz(i,j+1)<0.0) diffp=litho(k)%diffm
          diffm=litho(k)%diffa
          if(gz(i,j-1)<0.0) diffm=litho(k)%diffm
          halfyp=0.5*(diff*Flitho(i,j,k)+diffp*Flitho(i,j+1,k))
          halfym=0.5*(diff*Flitho(i,j,k)+diffm*Flitho(i,j-1,k))
          ax=halfxp*(gz(i+1,j)-gz(i,j))-halfxm*(gz(i,j)-gz(i-1,j))
          ay=halfyp*(gz(i,j+1)-gz(i,j))-halfym*(gz(i,j)-gz(i,j-1))
          sum=sum+ax/dx2+ay/dx2
        enddo
        ! In case of erosion limit the time step for stability
        if(sum<0.)then
          dt=min(dt,-topLay*0.75/sum)
        endif
        hrate(i,j)=sum
      enddo
    enddo
    dt=dt*0.1
    ! Adjust time step to match with display
    if(time+dt>display) dt=display-time
    ! Adjust time step to match with layer time
    if(time+dt>strati) dt=strati-time
    call mpi_allreduce(mpi_in_place,dt,1,mpi_double_precision,mpi_min,comm,ierr)

    ! Temporary elevations values

    tempz=0.
    do j=js(pet_id+1)+1,je(pet_id+1)-1
      do i=is(pet_id+1)+1,ie(pet_id+1)-1
        tempz(i,j)=gz(i,j)+hrate(i,j)*dt
      enddo
    enddo

    ! Update boundary
    call update_Bounds_z

    ! Ghost cells
    tempz(1,1)=2.*tempz(2,2)-tempz(3,3)
    tempz(1,ny)=2.*tempz(2,ny-1)-tempz(3,ny-2)
    tempz(nx,1)=2.*tempz(nx-1,2)-tempz(nx-2,3)
    tempz(nx,ny)=2.*tempz(nx-1,ny-1)-tempz(nx-2,ny-2)
    tempz(1,2:ny-1)=2.*tempz(2,2:ny-1)-tempz(3,2:ny-1)
    tempz(nx,2:ny-1)=2.*tempz(nx-1,2:ny-1)-tempz(nx-2,2:ny-1)
    tempz(2:nx-1,1)=2.*tempz(2:nx-1,2)-tempz(2:nx-1,3)
    tempz(2:nx-1,ny)=2.*tempz(2:nx-1,ny-1)-tempz(2:nx-1,ny-2)

    ! Solve sediment proportion equation
    tmpLitho=0.
    do k=1,nlith-1
      do j=js(pet_id+1)+1,je(pet_id+1)-1
        do i=is(pet_id+1)+1,ie(pet_id+1)-1
          diff=litho(k)%diffa
          if(tempz(i,j)<0.0) diff=litho(k)%diffm
          diffp=litho(k)%diffa
          if(tempz(i+1,j)<0.0) diffp=litho(k)%diffm
          diffm=litho(k)%diffa
          if(tempz(i-1,j)<0.0) diffm=litho(k)%diffm
          halfxp=0.
          if(tempz(i-1,j)>tempz(i+1,j))then
            halfxp=dt*(diff*Flitho(i,j,k)-diffm*Flitho(i-1,j,k))/(2.*dx2)
          else
            halfxp=dt*(diffp*Flitho(i+1,j,k)-diff*Flitho(i,j,k))/(2.*dx2)
          endif
          diffp=litho(k)%diffa
          if(tempz(i,j+1)<0.0) diffp=litho(k)%diffm
          diffm=litho(k)%diffa
          if(tempz(i,j-1)<0.0) diffm=litho(k)%diffm
          halfyp=0.
          if(tempz(i,j-1)>tempz(i,j+1))then
            halfyp=dt*(diff*Flitho(i,j,k)-diffm*Flitho(i,j-1,k))/(2.*dx2)
          else
            halfyp=dt*(diffp*Flitho(i,j+1,k)-diff*Flitho(i,j,k))/(2.*dx2)
          endif
          sum=halfxp*(tempz(i+1,j)-tempz(i-1,j))+halfyp*(tempz(i,j+1)-tempz(i,j-1))
          tmpLitho(i,j,k)=(sum+Flitho(i,j,k)*topLay)/(topLay+tempz(i,j)-gz(i,j))
        enddo
      enddo
    enddo

    ! Update boundary
    call update_Bounds_f

    ! Update elevation
    gz=tempz

    ! Get all lithology proportions
    do j=js(pet_id+1)+1,je(pet_id+1)-1
      do i=is(pet_id+1)+1,ie(pet_id+1)-1
        sum=0.
        do k=1,nlith-1
          sum=sum+tmpLitho(i,j,k)
          Flitho(i,j,k)=tmpLitho(i,j,k)
        enddo
        if(sum>=1.)then
          diff=1./sum
          do k=1,nlith-1
            Flitho(i,j,k)=diff*Flitho(i,j,k)
          enddo
          Flitho(i,j,nlith)=0.
        else
          Flitho(i,j,nlith)=1.-sum
        endif
      enddo
    enddo

    ! Update the layers and compute eroding cells new proportions
    do j=js(pet_id+1)+1,je(pet_id+1)-1
      do i=is(pet_id+1)+1,ie(pet_id+1)-1
        if(hrate(i,j)>0.)then
          do k=1,nlith
            nlitho(k)=Flitho(i,j,k)
          enddo
          th=hrate(i,j)*dt
        elseif(hrate(i,j)<=0.)then
          nlitho(1:nlith)=Flitho(i,j,1:nlith)
          th=dt*hrate(i,j)
        else
          nlitho(1:nlith)=0.
          th=0.
        endif
        call update_layers(i,j,nlitho,th)
        Flitho(i,j,1:nlith)=nlitho(1:nlith)
      enddo
    enddo

    ! Last we update the proportion for the ghost cells
    do k=1,nlith
      Flitho(1,1,k)=Flitho(2,2,k)
      Flitho(1,ny,k)=Flitho(2,ny-1,k)
      Flitho(nx,1,k)=Flitho(nx-1,2,k)
      Flitho(nx,ny,k)=Flitho(nx-1,ny-1,k)
      Flitho(1,2:ny-1,k)=Flitho(2,2:ny-1,k)
      Flitho(nx,2:ny-1,k)=Flitho(nx-1,2:ny-1,k)
      Flitho(2:nx-1,1,k)=Flitho(2:nx-1,2,k)
      Flitho(2:nx-1,ny,k)=Flitho(2:nx-1,ny-1,k)
    enddo

    return

  end subroutine solve_explicit_PDE
  !---------------------------------------------------------------------------
  subroutine update_layers(i,j,nlitho,th)

    integer::i,j,l,k,lid

    real(kind=8)::th,toplayh,layh,temph,sumh
    real(kind=8),dimension(nlith)::nlitho,addh,frac,top

    ! Change in top layer thickness
    toplayh=topLay+th

    ! Case 1: layer increases, we add lithology from top layer
    ! to the stratigraphic pile
    if(toplayh>topLay)then
      do k=1,nlith
        stratiLay(i,j,nlay,k)=stratiLay(i,j,nlay,k)+nlitho(k)*th
        nlitho(k)=nlitho(k)
      enddo
    endif

    ! Case 2: layer decreases, we adjust active layer composition
    ! and the stratigraphic pile
    if(toplayh<topLay)then
      temph=toplayh
      addh=0.
      ! Take sediment thickness as required
      layloop: do l=nlay,1,-1
        layh=0
        frac=0.
        do k=1,nlith
          layh=layh+stratiLay(i,j,l,k)
          frac(k)=stratiLay(i,j,l,k)
        enddo
        ! The layer exists (e.g. not a sedimentary hiatus)
        if(layh>0.)then
          frac=frac/layh
          ! Stratal layer too small
          if(layh<=topLay-temph)then
            do k=1,nlith
              addh(k)=addh(k)+stratiLay(i,j,l,k)
              temph=temph+stratiLay(i,j,l,k)
              stratiLay(i,j,l,k)=0.
            enddo
          ! Otherwise take a fraction of it
          else
            lid=l
            do k=1,nlith
              addh(k)=addh(k)+frac(k)*(topLay-temph)
              stratiLay(i,j,l,k)=stratiLay(i,j,l,k)-frac(k)*(topLay-temph)
            enddo
            temph=topLay
          endif
          if(temph>=topLay) exit layloop
        endif
      enddo layloop

      ! Update each lithology fraction
      sumh=0
      top=0
      do k=1,nlith
        top(k)=addh(k)+nlitho(k)*toplayh
        sumh=sumh+addh(k)+nlitho(k)*toplayh+stratiLay(i,j,lid,k)
        nlitho(k)=addh(k)+nlitho(k)*toplayh+stratiLay(i,j,lid,k)
      enddo
      nlitho=nlitho/sumh

      ! Update straigraphic top layer thickness
      do k=1,nlith
        stratiLay(i,j,lid,k)=top(k)+stratiLay(i,j,lid,k)-nlitho(k)*topLay
      enddo

    endif

    return

  end subroutine update_layers
  !---------------------------------------------------------------------------

end module diffusionSolve
