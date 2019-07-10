program EF_analysis

 implicit none

 complex(16), allocatable :: chi(:,:), chi_dt(:,:), chi_dR(:,:,:)
 real(8),     allocatable :: V(:,:,:), w(:), tdPES_BO(:), grid_value(:,:), G(:,:)
 integer                  :: n_dof, n_grid, nstates, time_steps

 integer                  :: i, j, nlines, io, time
 character(100)           :: pot_file, chi_file, nullo
 complex(16)              :: norm

 pot_file='EF_parameter_gV_1024'
 chi_file='EF_parameter_dpsi_1024'

 open(file=pot_file,unit=100,status='old',action='read')

   read(100,*) nullo, nstates
   read(100,*) nullo, n_dof
   read(100,*) nullo, n_grid

   allocate(V(nstates,nstates,n_grid), w(n_grid), grid_value(n_grid,n_dof), G(n_dof,n_dof))

   read(100,*)
   read(100,*)
   do i=1,n_dof
    read(100,*) j,G(:,i)
   enddo

   read(100,*)
   read(100,*)

   do i=1,n_grid
    read(100,*) nullo, j, w(i), grid_value(i,:), V(:,:,i)
   enddo

 close(100)

 open(file=chi_file,unit=100,status='old',action='read')

   nlines = 0
   do
    read(100,*,iostat=io)
    if(io/=0) exit
    nlines = nlines + 1
   enddo

 close(100)

 time_steps=(nlines-1)/(n_grid+1)

 allocate(chi(n_grid,nstates),chi_dt(n_grid,nstates),chi_dR(n_grid,n_dof,nstates))

 open(file=chi_file,unit=100,status='old',action='read')

 read(100,*)

  do time=1,time_steps
   do i=1,n_grid
    read(100,*) nullo, nullo, grid_value(i,:), chi(i,:), chi_dt(i,:), chi_dR(i,:,:)
   enddo
   read(100,*)
   do j=1,nstates
    norm=dot_product(chi(:,j),w(:)*chi(:,j))
    write(*,*) norm
   enddo
  enddo

 close(100)


end program EF_analysis
