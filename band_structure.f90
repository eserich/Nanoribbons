!#/cm/local/apps/gcc/7.2.0/bin/gfortran

program brick_BND


!! Gap of armchair nanoribbons is located at Gamma
!! The program calculates the energy eigenvalues of a brick model
!topologically equivalent to hexagonal nanoriboons and more general than
!ladder
!model

 implicit none


!Variables
 integer :: j , l
 integer :: number_rows , kpoints
 real :: onsite_term , hopping_term , energy , k
 real , dimension(2) :: kpath
 real , dimension(:) , allocatable :: k_vector



!Reading input
open( unit=4,File='input.dat',status='OLD')

read(4,'(3x,i1)') number_rows
read(4, '(8x,i2)') kpoints
read(4, '(6x, f4.1, f4.1)') kpath
read(4, '(2x , f4.1)') hopping_term
read(4, '(2x , f4.1)') onsite_term

close(4)


!Calculating k_vector

allocate(k_vector(1:kpoints+1))
call calculate_kpoints(kpoints, kpath, k_vector)



!Writing eigenvalues
open( unit=2, file='output.dat',status='NEW')

write(2,*) ' #TB brick model of Graphene nanoribbons'

do j=1,number_rows

  write(2,*) '# band'
  do l=1,kpoints+1
      call eigenvalues(j,number_rows,k_vector(l),onsite_term,hopping_term,energy)
      write(2,*) l-1,  -energy
  end do
  write(2,*)
  write(2,*) '# band'
  do l=1,kpoints+1
      call eigenvalues(j,number_rows,k_vector(l),onsite_term,hopping_term,energy)
      write(2,*) l-1,  energy
  end do
  write(2,*)

end do

close(2)







j=4
k=0

call eigenvalues(j,number_rows,k,onsite_term,hopping_term,energy)
print*, energy


contains

subroutine eigenvalues(j,number_rows,k,onsite_term,hopping_term,energy)

 integer, intent(in) :: j , number_rows
 real , intent(in) :: k
 real , intent(in) :: hopping_term, onsite_term
 real , intent(out) :: energy
 real , parameter :: pi=3.14159265


   energy=sqrt(onsite_term**2 + (hopping_term**2)*(1+4*(cos(pi*j/(number_rows+1)))**2+4*cos(k*pi)*cos(pi*j/(number_rows+1))))


end subroutine eigenvalues

subroutine calculate_kpoints(kpoints,kpath,k_vector)

 integer :: j
 integer , intent(in) :: kpoints
 real , dimension(2) , intent(in) :: kpath
 real , dimension(kpoints+1) , intent(out) :: k_vector

 do j=1,kpoints+1
   k_vector(j) = kpath(1) + (j-1)*(kpath(2)-kpath(1))/(kpoints)
 end do

end subroutine calculate_kpoints











end program brick_BND

