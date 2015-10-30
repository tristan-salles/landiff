!---------------------------------------------------------------------------
!             AUTHOR TRISTAN SALLES -- UNIVERSITY OF SYDNEY
!---------------------------------------------------------------------------
program multidiff

  use diffusionParam
  use diffusionSolve
  implicit none

  integer::opt,k,dsp
  logical::exist
  real(kind=8)::maxlith
  real(kind=8)::t1,t2

  ! MPI Initialization
  call mpi_init(ierr)
  comm = MPI_COMM_WORLD
  call mpi_comm_size(comm,npets,ierr)
  call mpi_comm_rank(comm,pet_id,ierr)
  t1=mpi_wtime()

  ! Read command line
  if(pet_id==0)then
    opt=iargc()
    if(opt<1)then
      print*,'Wrong command line you need to specify initial topographic file'
      return
    elseif(opt==1)then
      call getarg(1,Fnodes)
    endif

    ! Check input file
    inquire(file=Fnodes,exist=exist)
    if(.not.exist)then
      print*,'Looking for'
      print*,trim(Fnodes)
      print*,'Topographic file not found please specify the path again.'
      return
    endif
  endif
  call mpi_bcast(Fnodes,128,mpi_character,0,comm,ierr)

  ! Time definition
  nlith=3
  start=0.
  end=1000000.
  time=start
  dispstep=50000.
  stratistep=10000.

  ! Determine stratigraphic grid
  nbLay=int((end-start)/stratistep)+2
  bottomlay=3000.

  ! Determine grid size
  call get_grid_dimension

  ! Partition grid
  call grid_partition

  ! Define lithological parameters
  topLay=10.
  maxlith=0.
  if(allocated(litho))deallocate(litho)
  allocate(litho(nlith))
  do k=1,nlith
    litho(k)%compac=1.
    litho(k)%diffm=0.1*k !*0.5
    litho(k)%diffa=0.5*k
    maxlith=max(maxlith,litho(k)%diffa)
    maxlith=max(maxlith,litho(k)%diffm)
  enddo

  if(allocated(Flitho))deallocate(Flitho)
  allocate(Flitho(nx,ny,nlith))
  allocate(tmpLitho(nx,ny,nlith))
  Flitho(1:nx,1:ny,1)=0.3
  Flitho(1:nx,1:ny,2)=0.3
  Flitho(1:nx,1:ny,3)=0.4

  ! Create stratigraphic mesh

  if(allocated(stratiLay))deallocate(stratiLay)
  allocate(stratiLay(nx,ny,nbLay,nlith))
  if(allocated(stratLayID))deallocate(stratLayID)
  allocate(stratLayID(nx,ny))
  nlay=1
  stratiLay(1:nx,1:ny,1,1)=bottomlay*0.3-topLay
  stratiLay(1:nx,1:ny,1,2)=bottomlay*0.3-topLay
  stratiLay(1:nx,1:ny,1,3)=bottomlay*0.4-topLay
  nlay=2
  stratLayID(1:nx,1:ny)=1


  ! Time step for stability
  maxstep=real(0.25*int(dx2/(2.*maxlith)))
  dsp=0
  ! Display time
  call write_output(dsp)
  display=start+dispstep

  ! Stratigraphic layer time
  strati=start+stratistep

  do while(time<=end)
    ! Fully explicit method
    call solve_explicit_PDE

    ! Update time
    time=time+dt
    if(pet_id==0) print*,'time - ',time,dt
    ! Display if needed
    if(time==display)then
      dsp=dsp+1
      call write_output(dsp)
      display=display+dispstep
    endif

    ! Stratigraphy if needed
    if(time==strati)then
      nlay=nlay+1
      strati=strati+stratistep
    endif

  enddo

  t2=mpi_wtime()
  if(pet_id==0) print*,'Simulation done...',t2-t1

  call mpi_finalize(ierr)

  return

end program multidiff
