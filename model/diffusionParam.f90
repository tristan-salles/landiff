!---------------------------------------------------------------------------
!             AUTHOR TRISTAN SALLES -- UNIVERSITY OF SYDNEY
!---------------------------------------------------------------------------
module diffusionParam

  implicit none
  include 'mpif.h'

  ! MPI data
  logical::rowPart
  integer::comm,ierr,npets,pet_id
  integer,dimension(2)::NeighBor
  integer,dimension(:),allocatable::is,js,ie,je

  ! File
  character(len=128)::Fnodes

  ! Time parameters
  real(kind=8)::start,end,time
  real(kind=8)::display,dispstep
  real(kind=8)::strati,stratistep

  ! Grid parameters
  integer::nx,ny,xcell,ycell
  integer,dimension(:),allocatable::xpart,ypart
  integer::xdomains,ydomains
  real(kind=8)::xmin,ymin,xmax,ymax,dx,dx2
  real(kind=8),dimension(:,:),allocatable::gx,gy,gz,tempz
  real(kind=8),dimension(:,:,:),allocatable::tmpLitho


  ! Lithology
  real(kind=8)::topLay
  integer::nlith
  type lithology_parameters
     ! Compaction ratio
     real(kind=8)::compac
     ! Marine diffusivity
     real(kind=8)::diffm
     ! Aerial diffusivity
     real(kind=8)::diffa
  end type lithology_parameters
  type(lithology_parameters),dimension(:),allocatable::litho
  ! Fraction of lithology within the top layer
  real(kind=8),dimension(:,:,:),allocatable::Flitho

  ! Stratigraphic layer
  integer::nbLay,nlay
  real(kind=8)::bottomlay
  real(kind=8),dimension(:,:),allocatable::stratLayID
  real(kind=8),dimension(:,:,:,:),allocatable::stratiLay

  ! Time step
  real(kind=8)::maxstep,dt

contains

  !---------------------------------------------------------------------------
  subroutine get_grid_dimension

    integer::ios,n,nbnodes,i,j

    integer,parameter::maxrecs=1000000000

    real(kind=8)::xo,yo,h

    ! Read topographic file
    open(12,file=Fnodes,status="old",action="read",iostat=ios)

    nbnodes=0
    nodes_count: do n=1,maxrecs
      read(12,*,iostat=ios) xo,yo,h
      if(ios/=0)then
        xmax=xo
        ymax=yo
        exit nodes_count
      endif
      if(n==1)then
        xmin=xo
        ymin=yo
      elseif(n==2)then
        dx=xo-xmin
        dx2=dx*dx
      endif
      if(n==maxrecs)then
        print*,'Maximum number of nodes number exceeded in the regular input file'
        exit
      endif
      nbnodes=nbnodes+1
    enddo nodes_count
    rewind(12)

    ! Derived the grid dimensions
    nx=int((xmax-xmin)/dx)+1
    ny=int((ymax-ymin)/dx)+1
    ! Add the ghost cells
    nx=nx+2
    ny=ny+2

    ! Allocate nodes coordinates
    if(allocated(gx))deallocate(gx)
    if(allocated(gy))deallocate(gy)
    if(allocated(gz))deallocate(gz)
    allocate(gx(nx,ny),gy(nx,ny),gz(nx,ny),tempz(nx,ny))
    do j=2,ny-1
      do i=2,nx-1
        read(12,*)gx(i,j),gy(i,j),gz(i,j)
        gz(i,j)=gz(i,j) !+bottomlay
      enddo
    enddo

    ! Get ghost point values
    ! For x-coordinates
    gx(1,1:ny)=xmin-dx
    gx(nx,1:ny)=xmax+dx
    gx(1:nx,1)=gx(1:nx,2)
    gx(1:nx,ny)=gx(1:nx,ny-1)
    ! For y-coordinates
    gy(1:nx,1)=ymin-dx
    gy(1:nx,ny)=ymax+dx
    gy(1,1:ny)=gy(2,1:ny)
    gy(nx,1:ny)=gy(nx-1,1:ny)
    ! For z-coordinates
    gz(1,1)=2.*gz(2,2)-gz(3,3)
    gz(1,ny)=2.*gz(2,ny-1)-gz(3,ny-2)
    gz(nx,1)=2.*gz(nx-1,2)-gz(nx-2,3)
    gz(nx,ny)=2.*gz(nx-1,ny-1)-gz(nx-2,ny-2)
    gz(1,2:ny-1)=2.*gz(2,2:ny-1)-gz(3,2:ny-1)
    gz(nx,2:ny-1)=2.*gz(nx-1,2:ny-1)-gz(nx-2,2:ny-1)
    gz(2:nx-1,1)=2.*gz(2:nx-1,2)-gz(2:nx-1,3)
    gz(2:nx-1,ny)=2.*gz(2:nx-1,ny-1)-gz(2:nx-1,ny-2)
    xmax=xmax+dx
    xmin=xmin-dx
    ymax=ymax+dx
    ymin=ymin-dx

    ! Close topographic file
    close(12)

    return

  end subroutine get_grid_dimension
  !---------------------------------------------------------------------------
  subroutine grid_partition

    integer::part(npets)
    integer::k,remain,cum

    allocate(is(npets))
    allocate(ie(npets))
    allocate(js(npets))
    allocate(je(npets))
    is=0
    ie=0
    js=0
    je=0
    ! Column partitioning
    if(nx>=ny)then
      rowPart=.false.
      part=int(nx/npets)
      remain=nx-int(nx/npets)*npets
      do k=1,remain
        part(k)=part(k)+1
      enddo
      is(1)=1
      ie(npets)=nx
      cum=part(1)
      do k=2,npets
        ie(k-1)=cum+1
        is(k)=cum-1
        cum=cum+part(k)
      enddo
      js=1
      je=ny
    ! Row partitioning
    else
      rowPart=.true.
      part=int(ny/npets)
      remain=ny-int(ny/npets)*npets
      do k=1,remain
        part(k)=part(k)+1
      enddo
      js(1)=1
      je(npets)=ny
      cum=part(1)
      do k=2,npets
        je(k-1)=cum+1
        js(k)=cum-1
        cum=cum+part(k)
      enddo
      is=1
      ie=nx
    endif

    NeighBor(:) = MPI_Proc_null
    if(pet_id/=0) NeighBor(1)=pet_id-1
    if(pet_id/=npets-1) NeighBor(2)=pet_id+1

    return

  end subroutine grid_partition
  !---------------------------------------------------------------------------
  subroutine update_Bounds_z

    integer::flag
    integer,dimension(mpi_status_size)::status

    flag=1
    if(.not.rowPart)then
      ! Send my boundary to East and receive from West
      call mpi_sendrecv(tempz(ie(pet_id+1)-2,1:ny),ny,mpi_double_precision,neighBor(2),flag, &
        tempz(is(pet_id+1),1:ny),ny,mpi_double_precision,neighBor(1),flag,comm,status,ierr)
      ! Send my boundary to West and receive from East
      call mpi_sendrecv(tempz(is(pet_id+1)+2,1:ny),ny,mpi_double_precision,neighBor(1),flag, &
        tempz(ie(pet_id+1),1:ny),ny,mpi_double_precision,neighBor(2),flag,comm,status,ierr)
    endif

    flag = 2
    if(rowPart)then
      ! Send my boundary to North and receive from South
      call mpi_sendrecv(tempz(1:nx,je(pet_id+1)-2),nx,mpi_double_precision,neighBor(2),flag, &
        tempz(1:nx,js(pet_id+1)),nx,mpi_double_precision,neighBor(1),flag,comm,status,ierr)

      ! Send my boundary to South and receive from North
      call mpi_sendrecv(tempz(1:nx,js(pet_id+1)+2),nx,mpi_double_precision,neighBor(1),flag, &
        tempz(1:nx,je(pet_id+1)),nx,mpi_double_precision,neighBor(2),flag,comm,status,ierr)
    endif

    return

  end subroutine update_Bounds_z
  !---------------------------------------------------------------------------
  subroutine update_Bounds_f

    integer::flag,k
    integer,dimension(mpi_status_size)::status

    flag=1
    if(.not.rowPart)then
      do k=1,nlith-1
        ! Send my boundary to East and receive from West
        call mpi_sendrecv(tmpLitho(ie(pet_id+1)-2,1:ny,k),ny,mpi_double_precision,neighBor(2),flag, &
          tmpLitho(is(pet_id+1),1:ny,k),ny,mpi_double_precision,neighBor(1),flag,comm,status,ierr)
        ! Send my boundary to West and receive from East
        call mpi_sendrecv(tmpLitho(is(pet_id+1)+2,1:ny,k),ny,mpi_double_precision,neighBor(1),flag, &
          tmpLitho(ie(pet_id+1),1:ny,k),ny,mpi_double_precision,neighBor(2),flag,comm,status,ierr)
      enddo
    endif

    flag = 2
    if(rowPart)then
      do k=1,nlith-1
        ! Send my boundary to North and receive from South
        call mpi_sendrecv(tmpLitho(1:nx,je(pet_id+1)-2,k),nx,mpi_double_precision,neighBor(2),flag, &
          tmpLitho(1:nx,js(pet_id+1),k),nx,mpi_double_precision,neighBor(1),flag,comm,status,ierr)

        ! Send my boundary to South and receive from North
        call mpi_sendrecv(tmpLitho(1:nx,js(pet_id+1)+2,k),nx,mpi_double_precision,neighBor(1),flag, &
          tmpLitho(1:nx,je(pet_id+1),k),nx,mpi_double_precision,neighBor(2),flag,comm,status,ierr)
          if(pet_id==neighBor(2))print*,tmpLitho(1:nx,je(pet_id+1),k)
      enddo
    endif

    return

  end subroutine update_Bounds_f
  !---------------------------------------------------------------------------
  subroutine noblnk(string)

    integer::i,j,lg

    character(len=128)::string

    lg=len(string)
    do
       if(lg<=0 .or. string(lg:lg)/=' ') exit
       lg=lg-1
    enddo

    if(lg>0)then
       ! find first non-blank character
       i=1
       do
          if(i>lg .or. (string(i:i)/=' ' .and. string/=' ')) exit
          i=i+1
       enddo
       ! determine end of continuous (non-blank) string
       j=i
       do
          if(j>lg)then
             exit
          elseif(string(j:j)==' ')then
             exit
          elseif(string=='  ')then
             exit
          elseif(j==128)then
             exit
          endif
          j=j+1
       enddo
       ! j points to first blank position or past end of string; adjust to last
       ! non-blank position in string
       j=min(j-1,lg)
       string=string(i:j)
       if(j<len(string)) string(j+1:len(string))=' '
    else
       ! there were only blanks in string
       string=' '
    endif

    return

  end subroutine noblnk
  !---------------------------------------------------------------------------
  subroutine append_nb(stg1,i)

    integer::l1,l2,i

    character(len=128)::stg1,stg2

    l1=len_trim(stg1)
    write(stg2,'(i10)')i
    call noblnk(stg2)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2
    call noblnk(stg1)

    return

  end subroutine append_nb
  !---------------------------------------------------------------------------
  subroutine append_str(stg1,stg2)

    integer :: l1,l2
    character(len=128) :: stg1,stg2

    l1=len_trim(stg1)
    l2=len_trim(stg2)
    stg1(l1+1:l1+l2)=stg2
    call noblnk(stg1)

    return

  end subroutine append_str
  !---------------------------------------------------------------------------
  subroutine write_output(dsp)

    integer::dsp,i,j,ios,l1

    character(len=128)::Foutput,stg

    ! Create output file
    l1=len_trim(Fnodes)
    Foutput=Fnodes(1:l1-4)
    stg=''
    stg(1:2)='-p'
    call append_str(Foutput,stg)
    call append_nb(Foutput,pet_id)
    stg=''
    stg(1:1)='.'
    call append_str(Foutput,stg)
    call append_nb(Foutput,dsp)
    stg=''
    stg(1:4)='.csv'
    call append_str(Foutput,stg)

    ! Open output file
    open(12,file=Foutput,status="replace",action="write",iostat=ios)
    write(12,*)'x,y,z,elev,f1,f2,f3'
    do j=js(pet_id+1)+1,je(pet_id+1)-1
      do i=is(pet_id+1)+1,ie(pet_id+1)-1
        write(12,*)gx(i,j),',',gy(i,j),',',gz(i,j),',',gz(i,j),',',Flitho(i,j,1),',',Flitho(i,j,2),',',Flitho(i,j,3)
      enddo
    enddo
    close(12)

    return

  end subroutine write_output
  !---------------------------------------------------------------------------

end module diffusionParam
