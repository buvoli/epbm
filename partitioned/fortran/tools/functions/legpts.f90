function legpts(n,I_in) result(x)
! Returns legendre points on interval I by applying the recurrence relation for the Legendre polynomials and their derivatives to perform Newton iteration on the WKB approximation to the roots.
!This code was adapted from Chebfun (https://www.chebfun.org) legpts.m
! === Parameters ===========================================
! n - number of chebyshev points
! I_in - (array) 2x1 array containing inteval
! === Output ==============================================
! c - (array) nx1 array of chebyshev points


integer, intent(in) :: n
real(kind=8), intent(in), optional :: I_in(2)
real(kind=8) :: I(2), eps
real(kind=8), allocatable :: x(:), Pm1(:), Pm2(:), PPm1(:), PPm2(:), P(:), PP(:), dx(:), xp(:)
integer :: counter, j, npn, nnn

if(present(I_in)) then
    I = I_in
else
    I = (/ -1.d0, 1.d0 /)
endif

allocate(x(n))

if(n .eq. 1) then
  x(1) = 0.d0
  RETURN
endif


npn = ceiling(n / 2.d0) ! num positive nodes
nnn = floor(n / 2.0d0)  ! num negative nodes

allocate(xp(npn))
allocate(P(npn))
allocate(PP(npn))
allocate(Pm1(npn))
allocate(Pm2(npn))
allocate(PPm2(npn))
allocate(PPm1(npn))
allocate(dx(npn))

do j = 1, npn
    xp(j) = real(npn - (j - 1))
enddo

! form initial guess
xp = PI * (4.d0 * xp - 1.d0) / (4.d0 * n + 2.d0)
xp = (1.d0 - (n - 1.d0) / (8.d0 * n**3) - 1.d0 / (384.d0 * n**4) * (39.d0 - 28.d0 / sin(xp)**2) ) * cos(xp)

Pm2  = 1.d0
Pm1  = xp
PPm2 = 0.d0
PPm1 = 1.d0
eps  = .00000000000000023d0
dx   = 2.d0 * eps
counter = 0

! Newton Interations
do while ( maxval(abs(dx)) > eps .and. counter < 10)
    counter = counter + 1;
    do j = 1, n-1
        P    = ((2.d0 * j + 1.d0) * Pm1 * xp - j * Pm2) / (j + 1.d0);
        Pm2  = Pm1;
        Pm1  = P;
        PP   = ((2.d0 * j + 1.d0) * (Pm2 + xp * PPm1) - j * PPm2) / (j + 1.d0);
        PPm2 = PPm1;
        PPm1 = PP;
    enddo
    ! Newton correction and update
    dx = -P / PP;
    xp = xp + dx;
    ! reinitialize
    Pm2  = 1;
    Pm1  = xp;
    PPm2 = 0;
    PPm1 = 1;
enddo

! reflect positive points
x(n-npn+1:n) = xp
do j = 1, nnn
    x(j) = -1 * xp(npn - j + 1)
end do

! scale according to interval
x = (I(2) - I(1))/2.d0 * x + (I(2) + I(1))/2.d0

end function legpts
