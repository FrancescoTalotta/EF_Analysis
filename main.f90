program EF_analysis

 implicit none

 complex(16), allocatable :: chi(:,:), chi_dt(:,:), chi_dR(:,:,:), chi_dtdR(:,:,:), chi2(:)
 real(8),     allocatable :: V(:,:,:), w(:), tdPES_BO(:), grid_value(:,:), G(:,:)
 integer,     allocatable :: n_gridbasis(:)
 real(8)                  :: Tmax, deltaT
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
  
   allocate(n_gridbasis(n_dof))

   read(100,*) nullo, n_gridbasis(:)
   read(100,*) nullo, nullo, Tmax, deltaT 

   allocate(V(nstates,nstates,n_grid), w(n_grid), grid_value(n_grid,n_dof), G(n_dof,n_dof))

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

 time_steps=nint(Tmax/deltaT)+1

 allocate(chi(n_grid,nstates),chi_dt(n_grid,nstates),chi_dR(n_grid,n_dof,nstates),chi_dtdR(n_grid,n_dof,nstates),chi2(n_grid))

 open(file=chi_file,unit=100,status='old',action='read')

 read(100,*)

  do time=1,time_steps

     do i=1,n_grid
      read(100,*) nullo, nullo, grid_value(i,:), chi(i,:), chi_dt(i,:), chi_dR(i,:,:), chi_dtdR(i,:,:)
     enddo

     read(100,*)

     chi2=(0.d0,0.d0)  
     do j=1,nstates
      chi2(:)=chi2(:)+conjg(chi(:,j))*chi(:,j)
     enddo
     
     call tdpes() 

  enddo

 close(100)

end program EF_analysis


!subroutine chi2(chi,n_grid,nstates,n_dof,)
!
!  implicit none
!
!  complex(16), intent(in) :: chi
!  integer,     intent(in) :: n_grid, nstates, n_dof
!
!end subroutine chi2

!subroutine phase(chi, chi_dR, n_grid, nstates, n_dof, S, S_dR)
!
!  implicit none 
!
!      integer, intent(in) :: n_grid, nstates, n_dof
!  complex(16), intent(in) :: chi(n_grid,nstates), chi_dR(n_grid,n_dof,nstates) 
!     real(8), intent(out) :: S(n_grid), S_dR(n_grid) 
!    
!                  integer :: i
!
!  do i=1,n_grid 
!  enddo
!
!end subroutine phase
