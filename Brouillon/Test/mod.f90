module mesh
implicit none
type :: geometry
real(8) :: x0, x1, dx ! coordinates of origin and grid size
integer :: nx ! number of grid points
real(8), dimension(:), pointer :: xgrid ! coordinates of points
end type geometry
contains

subroutine create(geom,x0,nx,dx)
!f2py integer, intent(out) :: geom
type(geometry), pointer :: geom
real(8), intent(in) :: x0, dx
integer, intent(in) :: nx
integer :: i
allocate(geom)
geom%x0=x0; geom%x1=x0+nx*dx; geom%dx=dx; geom%nx=nx
allocate(geom%xgrid(nx))
do i=1,nx
geom%xgrid(i)=geom%x0+(i-1)*geom%dx
end do
end subroutine create

subroutine view
!f2py integer, intent(in) :: this
type(geometry),pointer :: this
!print*, ’nx = ’, this%nx
!print*, this%x0,this%x1
!print*, this%xgrid(:)
end subroutine view
end module mesh
