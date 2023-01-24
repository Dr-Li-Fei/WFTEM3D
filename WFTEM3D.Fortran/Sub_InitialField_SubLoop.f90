    !+ Calculate the subloop initial field.
    !
    Subroutine Sub_InitialField_SubLoop
    !
    ! In: i0,j0,k0,XI,YJ,ZK,model_EC,dx,dy,dz,subloop_x,subloop_y,L_subloop,
    !     Alpha (NOTE: model_EC is not necessary)
    ! Out: EX_subloop,EY_subloop,EZ_subloop,HX_subloop,HY_subloop,HZ_subloop,
    !      t0_E,t0_H
    ! Description:
    ! Calculate the subloop initial field.
    !
    ! Method:
    ! Calculate the whole-space initial field excited by magnetic dipole 
    ! sources at an initial time after the current is switched off. 
    ! See 3D Finite-difference Transient Electromagnetic Modeling with 
    ! Whole-space Initial Field for detail.
    !
    ! History:
    ! Version    Date    Comment
    ! -------    ----    -------
    ! 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
    !
    ! Code Description:
    ! Language: Fortran 90.
    ! Software Standards: "European Standards for Writing and
    ! Documenting Exchangeable Fortran 90 Code".
    !
    ! Declarations:
    ! Modules used:
    use Mod_Parameters
    ! Imported Parameters/Variables with intent:
    ! Permeability_vac     ! Magnetic permeability of vacuum.
    ! i0, j0, k0           ! Tx loop center is at the bottom face center
                           ! of the cell (i0,j0,k0).
    ! XI, YJ, ZK           ! Number of cells in the x-, y-, and z-directions.
    ! model_EC             ! Model conductivity.
    ! dx, dy, dz           ! Grid size in the x-, y-, and z-directions.
    ! subloop_x,subloop_y  ! Coordinate of subloop center.
    ! L_subloop            ! Length of subloops.
    ! Alpha                ! Coefficient in the time step size calculation equation.
    ! EX_subloop(:,:,:,:)  ! Subloop initial field Ex.
    ! EY_subloop(:,:,:,:)  ! Subloop initial field Ey.
    ! EZ_subloop(:,:,:,:)  ! Subloop initial field Ez.
    ! HX_subloop(:,:,:,:)  ! Subloop initial field Hx.
    ! HY_subloop(:,:,:,:)  ! Subloop initial field Hy.
    ! HZ_subloop(:,:,:,:)  ! Subloop initial field Hz.   
    ! t0_E                 ! Initial time of electric field.
    ! t0_H                 ! Initial time of magnetic field.
    ! real(8)::dt(0:1)     ! Time step size.
    implicit none
    ! Subroutine arguments:
    integer::i,j,k           ! Temporary loop variables.
    real(8)::coord_x,coord_y,coord_z ! Coordinates of the node.
    real(8)::coord_r         ! Distance between the grid node and Tx center.
    real(8)::u               ! Temporary variable used to calculate initial field.
    real(8)::initialfield_EC ! The conductivity used to calculate initial field.
    real(8)::Txmm            ! Tx magnetic moment.
    !- End of header --------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [1.0] Set the conductivity used to calculate initial field:
    !-----------------------------------------------------------------------------
    ! NOTE: The initial field is independent of the model, and can be calculated
    !       using an artificial high conductivity (e.g. initialfield_EC = 100),
    !       here the conductivity of cell (i0,j0,k0+1) is used by default:
    initialfield_EC = model_EC(i0-1,j0-1,k0)
    ! Tx magnetic moment (current is 1 A by default):
    Txmm = L_subloop*L_subloop
    !-----------------------------------------------------------------------------
    ! [2.0] Calculate initial time:
    !-----------------------------------------------------------------------------   
    ! Calculate the initial time of electric field using equation 26 in paper:
    t0_E = 1.13D0*Permeability_vac*model_EC(I0-1,J0-1,K0-1)*DZ(K0-1)**2
    ! Calculate the time step size using equation 28 in paper:
    dt(0) = Alpha*DZ(K0-1)*sqrt(Permeability_vac*model_EC(I0-1,J0-1,K0-1)*t0_E/6.0D0)
    ! Calculate the initial time of magnetic field using equation 27 in paper:
    t0_H = t0_E+dt(0)/2.0d0
    !-----------------------------------------------------------------------------
    ! [3.0] Calculate subloop initial field using equations 25 in paper:
    !-----------------------------------------------------------------------------    
    ! Calculate subloop initial field Ex:
    do k = 1,ZK-1
        do j = 1,YJ-1
            do i = 0,XI-1
                coord_x = coor_x_calc(i,1)+subloop_x
                coord_y = coor_y_calc(j,1)+subloop_y
                coord_z = coor_z_calc(k,1)
                coord_r = sqrt(coord_x**2+coord_y**2+coord_z**2)
                u = coord_r*sqrt(2.0d0*pi*initialfield_EC/(1.0d7*t0_E))
                EX_subloop(i,j,k,0) = sqrt(2.0d0/pi)*Txmm*coord_y*u**5 &
                /(4.0d0*pi*coord_r**5*initialfield_EC*exp(u**2/2.0d0))    
            end do
        end do
    end do
    ! Calculate subloop initial field Ey:
    do k = 1,ZK-1
        do j = 0,YJ-1
            do i = 1,XI-1
                coord_x = coor_x_calc(i,2)+subloop_x
                coord_y = coor_y_calc(j,2)+subloop_y
                coord_z = coor_z_calc(k,2)
                coord_r = sqrt(coord_x**2+coord_y**2+coord_z**2)
                u = coord_r*sqrt(2.0d0*pi*initialfield_EC/(1.0d7*t0_E))             
                EY_subloop(i,j,k,0) = sqrt(2.0d0/pi)*Txmm*coord_x*u**5 &
                /(4.0d0*pi*coord_r**5*initialfield_EC*exp(u**2/2.0d0))
            end do
        end do
    end do
    ! Calculate subloop initial field Ez:    
    EZ_subloop(1:XI-1,1:YJ-1,0:ZK-1,0) = 0.0d0
    ! Calculate subloop initial field Hx:
    do k = 0,ZK-1
        do j = 0,YJ-1  
            do i = 0,XI
                coord_x = coor_x_calc(i,4)+subloop_x
                coord_y = coor_y_calc(j,4)+subloop_y
                coord_z = coor_z_calc(k,4)
                coord_r = sqrt(coord_x**2+coord_y**2+coord_z**2)
                u = coord_r*sqrt(2.0d0*pi*initialfield_EC/(1.0d7*t0_H))
                HX_subloop(i,j,k,0) = Txmm*coord_z*coord_x             &
                /(4.0d0*pi*coord_r**5)*(3.0d0*erf(u)-sqrt(2.0d0/pi)    &
                *u*(3.0d0+U**2)/exp(u**2/2.0d0))
            end do
        end do
    end do
    ! Calculate subloop initial field Hy:
    do k = 0,ZK-1
        do j = 0,YJ
            do i = 0,XI-1      
                coord_x = coor_x_calc(i,5)+subloop_x
                coord_y = coor_y_calc(j,5)+subloop_y
                coord_z = coor_z_calc(k,5)
                coord_r = sqrt(coord_x**2+coord_y**2+coord_z**2)
                u = coord_r*sqrt(2.0d0*pi*initialfield_EC/(1.0d7*t0_H))
                HY_subloop(i,j,k,0) = Txmm*coord_z*coord_y             &
                /(4.0d0*pi*coord_r**5)*(3.0d0*erf(u)-sqrt(2.0d0/pi)    &
                *u*(3.0d0+U**2)/exp(u**2/2.0d0))
            end do
        end do
    end do
    ! Calculate subloop initial field Hz:
    do k = 0,ZK
        do j = 0,YJ-1
            do i = 0,XI-1
                coord_x = coor_x_calc(i,6)+subloop_x
                coord_y = coor_y_calc(j,6)+subloop_y
                coord_z = coor_z_calc(k,6)
                coord_r = sqrt(coord_x**2+coord_y**2+coord_z**2)
                u = coord_r*sqrt(2.0d0*pi*initialfield_EC/(1.0d7*t0_H))
                if (coord_r < 0.001*dz(k0)) then
                   coord_r = 0.001*dz(k0) ! Avoid singularity.
                end if
                HZ_subloop(i,j,k,0) = Txmm/(4.0d0*pi*coord_r**5)             &
                *((2.0d0*coord_z**2-sqrt(coord_x**2+coord_y**2))*erf(u)      &
                -(2.0d0*coord_z**2-sqrt(coord_x**2+coord_y**2)*(1.0d0+u**2)) &
                *sqrt(2.0d0/pi)/exp(u**2/2.0d0))
            end do
        end do
    end do

    contains
    !-----------------------------------------------------------------------------
    ! [4.0] SubFunctions used to calculate coordinates in equations 25 in paper:
    !-----------------------------------------------------------------------------
    ! Calculate the value of x coordinate:
    Function coor_x_calc(x,m)
    integer:: x,m,i
    real(8):: coor_x_calc
    coor_x_calc = 0.0
    if (m == 2.or.m == 3.or.m == 4)then        
        if (x <= i0-1) then
            do i = x,i0-1
                coor_x_calc = coor_x_calc+dx(i)
            end do
            coor_x_calc = coor_x_calc-5.0d-1*dx(i0-1)
            coor_x_calc = -coor_x_calc
        elseif(x > i0-1) then
            do i = i0-1,x-1
                coor_x_calc = coor_x_calc+dx(i)
            end do
            coor_x_calc = coor_x_calc-5.0d-1*dx(i0-1)
        end if
    else
        if (x <= i0-1) then
            do i = x,i0-1
                coor_x_calc = coor_x_calc+dx(i)
            end do
            coor_x_calc = coor_x_calc-5.0d-1*(dx(x)+dx(i0-1))
            coor_x_calc = -coor_x_calc
        elseif(x > i0-1) then
            do i = i0-1,x                                                      
                coor_x_calc = coor_x_calc+dx(i)
            end do
            coor_x_calc = coor_x_calc-5.0d-1*(dx(i0-1)+dx(x))
        end if
    end if
    End Function
    ! Calculate the value of y coordinate:
    Function coor_y_calc(y,m)
    integer:: y,m,i
    real(8):: coor_y_calc
    coor_y_calc = 0.0
    if (m == 1.or.m == 3.or.m == 5)then
        if (y <= j0-1) then
            do i = y,j0-1
                coor_y_calc = coor_y_calc+dy(i)
            end do
            coor_y_calc = coor_y_calc-5.0d-1*dy(j0-1)
            coor_y_calc = -coor_y_calc
        elseif(y > j0-1) then
            do i = j0-1,y-1
                coor_y_calc = coor_y_calc+dy(i)
            end do
            coor_y_calc = coor_y_calc-5.0d-1*dy(j0-1)
        end if
    else
        if (y <= j0-1) then
            do i = y,j0-1
                coor_y_calc = coor_y_calc+dy(i)
            end do
            coor_y_calc = coor_y_calc-5.0d-1*(dy(y)+dy(j0-1))
            coor_y_calc = -coor_y_calc
        elseif(y > j0-1) then
            do i = j0-1,y
                coor_y_calc = coor_y_calc+dy(i)
            end do
            coor_y_calc = coor_y_calc-5.0d-1*(dy(j0-1)+dy(y))
        end if
    end if
    End Function
    ! Calculate the value of z coordinate:
    Function coor_z_calc(z,m)
    integer:: z,m,i
    real(8):: coor_z_calc
    coor_z_calc = 0.0
    if (m == 3.or.m == 4.or.m == 5)then
        if (z <= k0-1) then
            do i = z,k0-1
                coor_z_calc = coor_z_calc+dz(i)
            end do
            coor_z_calc = coor_z_calc-0.5*dz(z)
            coor_z_calc = -coor_z_calc
        elseif(z > k0-1) then
            do i = k0,z
                coor_z_calc = coor_z_calc+dz(i)
            end do
            coor_z_calc = coor_z_calc-0.5*dz(z)
        end if
    else
        if (z <= k0-1) then
            do i = z,k0-1
                coor_z_calc = coor_z_calc+dz(i)
            end do
            coor_z_calc = -coor_z_calc
        elseif(z > k0-1) then
            do i = k0,z-1
                coor_z_calc = coor_z_calc+dz(i)
            end do
        end if
    end if
    End Function
    EndSubroutine Sub_InitialField_SubLoop   
