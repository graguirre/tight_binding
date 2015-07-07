! =================================
!    Calculate all the distances
! =================================
SUBROUTINE dists(nx,ny,XYZ1,XYZ2,vec_dist)
implicit none
 INTEGER:: nx, ny, cont
 REAL*8,DIMENSION(nx,ny),intent(in):: XYZ1,XYZ2
 REAL*8,DIMENSION(nx*nx),intent(out):: vec_dist
 REAL*8,DIMENSION(3):: posi, posj, aux
 INTEGER::i,j
! REAL*8,DIMENSION(3):: zn

!zn = (/ 0.d0, 0.d0, 1.d0 /)
cont = 0
!$omp parallel default (shared) private (i)
!$omp do
do i = 1,nx
    posi = XYZ1(i,:)
    do j = 1,nx
      cont = cont + 1
      posj = XYZ2(j,:)
      aux = posi-posj
!      aux = aux - ( dot_product(aux,zn) * zn )
      vec_dist(cont) = sqrt(dot_product(aux, aux))
    end do
end do
!$omp end do
!$omp end parallel
END SUBROUTINE dists


! =========================
!   Numbering neighbours
! =========================
SUBROUTINE vecin(nx,ny,n,XYZ1,XYZ2,vec_dist,V,eps)
implicit none
 INTEGER:: nx, ny, n
 REAL*8,DIMENSION(nx,ny),intent(in):: XYZ1, XYZ2
 REAL*8,DIMENSION(n),intent(in):: vec_dist
 REAL*8,intent(in):: eps
 INTEGER,DIMENSION(nx,nx),intent(out):: V
 REAL*8,DIMENSION(3):: posi, posj, aux
 REAL*8:: r
 INTEGER:: i, j, k

do i = 1,nx
  do j = 1,nx
    V(i,j) = 0
  end do !j
end do !i

!$omp parallel default (shared) private (i)
!$omp do
do i = 1,nx
  posi = XYZ1(i,:)
  do j = 1,nx
    posj = XYZ2(j,:)
    aux = posi-posj
    r = sqrt(dot_product(aux, aux))
    do k = 1,n
      if ( r > vec_dist(k)-eps .and. r < vec_dist(k)+eps ) then
        V(i,j)=k-1
      end if
    end do !k
  end do !j
end do !i
!$omp end do
!$omp end parallel
END SUBROUTINE vecin



! =================================================
!    Computes the mean value of a given operator
!                    <vec| M |vec>
! =================================================
SUBROUTINE mean_val(nx, ny, vec, Mat, m)
implicit none
INTEGER:: nx, ny
COMPLEX*16:: vec(1, nx)
COMPLEX*16:: aux(1, nx)
COMPLEX*16:: Mat(nx, ny)
COMPLEX*16:: A(1,1)
COMPLEX*16, intent(out):: m

if (nx .ne. ny) then
  print*, 'FORTRAN: Operator is not a square matrix'
  stop
end if

aux = matmul(vec, mat)
A = matmul(aux,transpose(conjg(vec)))
m = A(1,1)
if (aimag(m) > 0.00001) then
  print*, 'FORTRAN: Mean val is complex'
  stop
end if
END SUBROUTINE mean_val


! =========
!    DOS
! =========
SUBROUTINE densityofstates(ne,nek,Ek,E,d,Y)
implicit none
 ! E: vector of energies in which we want to calculate the DOS
 ! Ek: collection of all eigenvals for every k
 ! Y: vector of DOS for each E
 INTEGER,intent(in):: ne, nek
 REAL*8,intent(in):: d
 REAL*8,DIMENSION(ne),intent(in):: E
 REAL*8,DIMENSION(nek),intent(in):: Ek
 REAL*8,DIMENSION(ne),intent(out):: Y
 ! dummy
 INTEGER:: i,j
 REAL*8:: sumando,de

do i=1,ne
  Y(i) = 0.d0
  do j=1,nek
    de = E(i)-Ek(j)
    if (de < 0.5) then
      sumando = d/((de)**2 + d*d)
      Y(i) = Y(i) + sumando
    end if
  end do
end do
END SUBROUTINE densityofstates
