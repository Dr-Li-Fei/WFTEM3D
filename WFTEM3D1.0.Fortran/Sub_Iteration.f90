    !+ Iteration and output results.
    !
    Subroutine Sub_Iteration
    !
    !-----------------------------------------------------------------------------
    ! In: i0,j0,k0,XI,YJ,ZK,dx,dy,dz,Rx_st,Rx_end,iter_n_max,Alpha,model_EC,
    !     EX,EY,EZ,HX,HY,HZ,t0_E,t0_H
    ! Out: t_iteration_H,DBZ_Rx
    ! Description:
    ! Iteration and output results.
    !
    ! Method:
    ! Step Maxwell¡¯s equations in time using a staggered grid and 
    ! a modified DuFort-Frankel method. See 3D Finite-difference Transient 
    ! Electromagnetic Modeling with Whole-space Initial Field for detail.
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
    ! dx, dy, dz           ! Grid size in the x-, y-, and z-directions.
    ! Rx_st, Rx_end        ! Receivers are at the bottom face center of cells
                           ! from (i0,Rx_st,k0) to (i0,Rx_end,k0).
    ! iter_n_max           ! Maximum number of iterations.
    ! Alpha                ! Coefficient in the time step size calculation equation.
    ! model_EC             ! Model conductivity.
    ! EX(:,:,:,:)          ! Ex.
    ! EY(:,:,:,:)          ! Ey.
    ! EZ(:,:,:,:)          ! Ez.
    ! HX(:,:,:,:)          ! Hx.
    ! HY(:,:,:,:)          ! Hy.
    ! HZ(:,:,:,:)          ! Hz.
    ! t0_E                 ! Initial time of electric field.
    ! t0_H                 ! Initial time of magnetic field.
    ! DBZ_Rx(:,:)          ! dBz/dt at receivers.
    ! t_iteration_E(:)     ! Time of electric field.
    ! t_iteration_H(:)     ! Time of magnetic field.
    ! EC_x(:,:,:)          ! Conductivity of nodes.
    ! EC_y(:,:,:)          ! ditto.
    ! EC_z(:,:,:)          ! ditto.
    ! dt(0:1)              ! Time step size.
    implicit none
    ! Subroutine arguments:
    integer::i,j,k,p                  ! Temporary loop variables.
    integer::iter_n                   ! Temporary iteration variable.
    character(len = 256)::n_Rx        ! Number of receivers.
    real(8)::gamma                    ! Artificial dielectric constant.
    !- End of header --------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [1.0] Calculate conductivity of grid nodes:
    !-----------------------------------------------------------------------------     
    ! Equation 13 in paper:
    do k = 1,ZK-1
        do j = 1,YJ-1 
            do i = 0,XI-1
                EC_x(i,j,k) = (dy(j-1)*dz(k-1)*model_EC(i,j-1,k-1) &
                +dy(j)*dz(k-1)*model_EC(i,j,k-1)                   &
                +dy(j-1)*dz(k)*model_EC(i,j-1,k)                   &
                +dy(j)*dz(k)*model_EC(i,j,k))                      &
                /((dy(j-1)+dy(j))*(dz(k-1)+dz(k)))
            end do
        end do
    end do
    ! Equation 15 in paper:    
    do k = 1,ZK-1
        do j = 0,YJ-1 
            do i = 1,XI-1    
                EC_y(i,j,k) = (dx(i-1)*dz(k-1)*model_EC(i-1,j,k-1) &
                +dx(i)*dz(k-1)*model_EC(i,j,k-1)                   &
                +dx(i-1)*dz(k)*model_EC(i-1,j,k)                   &
                +dx(i)*dz(k)*model_EC(i,j,k))                      &
                /((dx(i-1)+dx(i))*(dz(k-1)+dz(k)))
            end do 
        end do
    end do
    ! Equation 17 in paper:
    do k = 0,ZK-1
        do j = 1,YJ-1 
            do i = 1,XI-1
                EC_z(i,j,k) = (dx(i-1)*dy(j-1)*model_EC(i-1,j-1,k) &
                +dx(i)*dy(j-1)*model_EC(i,j-1,k)                   &
                +dx(i-1)*dy(j)*model_EC(i-1,j,k)                   &
                +dx(i)*dy(j)*model_EC(i,j,k))                      &
                /((dx(i-1)+dx(i))*(dy(j-1)+dy(j)))
            end do
        end do
    end do
    !-----------------------------------------------------------------------------
    ! [2.0] Get initial time:
    !-----------------------------------------------------------------------------
    t_iteration_E(0) = t0_E !Initial time of electric field.
    t_iteration_H(0) = t0_H !Initial time of magnetic field.
    ! Calculate dBz/dt at receivers at the initial time using equation 21 in paper:
    do p = Rx_st,Rx_end
        DBZ_Rx(0,p-Rx_st) = ((EX(i0-1,p,k0,0)-EX(i0-1,p-1,k0,0))/dy(p-1) &
        -(EY(i0,p-1,k0,0)-EY(i0-1,p-1,k0,0))/dx(i0-1))
    end do
    !-----------------------------------------------------------------------------
    ! [3.0] Iteration:
    !-----------------------------------------------------------------------------
    iter_n = 1
    do while (iter_n <= iter_n_max-1)
        print*, iter_n
        ! Calculate the electric field time of the iter_nth iteration:
        t_iteration_E(iter_n) = t_iteration_E(iter_n-1)+dt(0)
        ! Calculate the time step size using equation 23 in paper:       
        dt(1) = Alpha*dz(k0-1)*sqrt(Permeability_vac* &
        model_EC(i0-1,j0-1,k0-1)*t_iteration_E(iter_n)/6.0D0)
        ! Calculate the magnetic field time of the iter_nth iteration:
        t_iteration_H(iter_n) = t_iteration_H(iter_n-1)+(dt(0)+dt(1))/2.0d0
        ! Calculate artificial dielectric constant using equation 22 in paper:
        gamma = (4.0d0/Permeability_vac)*(dt(0)/dz(K0-1))**2
        ! Update the value of Ex using equation 12 in paper:
        do k = 1,ZK-1
            do j = 1,YJ-1 
                do i = 0,XI-1
                    EX(i,j,k,1) = (2.0d0*gamma-EC_x(i,j,k)*dt(0)) &
                    /(2.0d0*gamma+EC_x(i,j,k)*dt(0))*EX(i,j,k,0)  &                          
                    +4.0d0*dt(0)/(2.0d0*gamma+EC_x(i,j,k)*dt(0))  & 
                    *((HZ(i,j,k,0)-HZ(i,j-1,k,0))/(dy(j-1)+dy(j)) &                          
                    -(HY(i,j,k,0)-HY(i,j,k-1,0))/(dz(k-1)+dz(k))) 
                end do
            end do	   
        end do
        ! Update the value of Ey using equation 14 in paper:
        do k = 1,ZK-1
            do j = 0,YJ-1 
                do i = 1,XI-1             
                    EY(i,j,k,1) = (2.0d0*gamma-EC_y(i,j,k)*dt(0)) &
                    /(2.0d0*gamma+EC_y(i,j,k)*dt(0))*EY(i,j,k,0)  &
                    +4.0d0*dt(0)/(2.0d0*gamma+EC_y(i,j,k)*dt(0))  &
                    *((HX(i,j,k,0)-HX(i,j,k-1,0))/(dz(k-1)+dz(k)) &            
                    -(HZ(i,j,k,0)-HZ(i-1,j,k,0))/(dx(i-1)+dx(i)))	         
                end do
            end do	   
        end do
        ! Update the value of Ez using equation 16 in paper:
        do k = 0,ZK-1     
            do j = 1,YJ-1      
                do i = 1,XI-1             
                    EZ(i,j,k,1) = (2.0d0*gamma-EC_z(i,j,k)*dt(0)) &
                    /(2.0d0*gamma+EC_z(i,j,k)*dt(0))*EZ(i,j,k,0)  &
                    +4.0d0*dt(0)/(2.0d0*gamma+EC_z(i,j,k)*dt(0))  &
                    *((HY(i,j,k,0)-HY(i-1,j,k,0))/(dx(i-1)+dx(i)) &       
                    -(HX(i,j,k,0)-HX(i,j-1,k,0))/(dy(j-1)+dy(j)))

                end do
            end do
        end do 
        ! Dirichlet boundary condition:
        do i = 0,XI-1
            do j = 0,YJ
                EX(i,j,0,1) = 0.0d0
                EX(i,j,ZK,1) = 0.0d0
            end do
        end do
        do i = 0,XI-1
            do k = 0,ZK
                EX(i,0,k,1) = 0.0d0
                EX(i,YJ,k,1) = 0.0d0
            end do
        end do
        do i = 0,XI
            do j = 0,YJ-1
                EY(i,j,0,1) = 0.0d0
                EY(i,j,ZK,1) = 0.0d0
            end do
        end do
        do k = 0,ZK
            do j = 0,YJ-1
                EY(0,j,k,1) = 0.0d0
                EY(XI,j,k,1) = 0.0d0
            end do
        end do
        do i = 0,XI
            do k = 0,ZK-1
                EZ(i,0,k,1) = 0.0d0
                EZ(i,YJ,k,1) = 0.0d0
            end do
        end do
        do k = 0,ZK-1
            do j = 0,YJ
                EZ(0,j,k,1) = 0.0d0
                EZ(XI,j,k,1) = 0.0d0
            end do
        end do  
        ! Update the value of Hx using equation 18 in paper:
        do k = 0,ZK-1
            do j = 0,YJ-1
                do i = 0,XI
                    HX(i,j,k,1) = HX(i,j,k,0)           &
                    +(dt(0)+dt(1))/(2*Permeability_vac) &
                    *((EY(i,j,k+1,1)-EY(i,j,k,1))/dz(k) &
                    -(EZ(i,j+1,k,1)-EZ(i,j,k,1))/dy(j))
                end do
            end do	   
        end do
        ! Update the value of Hy using equation 19 in paper:
        do k = 0,ZK-1
            do j = 0,YJ
                do i = 0,XI-1
                    HY(i,j,k,1) = HY(i,j,k,0)           &
                    +(dt(0)+dt(1))/(2*Permeability_vac) &		
                    *((EZ(i+1,j,k,1)-EZ(i,j,k,1))/dx(i) &
                    -(EX(i,j,k+1,1)-EX(i,j,k,1))/dz(k))
                end do
            end do	   
        end do
        ! Update the value of Hz using equation 20 in paper, start from 
        ! the bottom boundary and the top boundary of the model where Hz=0:
        HZ(0:XI-1,0:YJ-1,0,1) = 0
        HZ(0:XI-1,0:YJ-1,ZK,1) = 0  
        do k = ZK,k0+1,-1
            do j = 0,YJ-1        
                do i = 0,XI-1
                    HZ(i,j,k-1,1) = HZ(i,j,k,1)                    &
                    +dz(k-1)/dx(i)*(HX(i+1,j,k-1,1)-HX(i,j,k-1,1)) &
                    +dz(k-1)/dy(j)*(HY(i,j+1,k-1,1)-HY(i,j,k-1,1))
                end do
            end do
        end do
        do k = 1,k0-1
            do j = 0,YJ-1        
                do i = 0,XI-1
                    HZ(i,j,k,1) = HZ(i,j,k-1,1)                    &
                    -dz(k-1)/dx(i)*(HX(i+1,j,k-1,1)-HX(i,j,k-1,1)) &
                    -dz(k-1)/dy(j)*(HY(i,j+1,k-1,1)-HY(i,j,k-1,1))
                end do
            end do
        end do
        ! Calculate dBz/dt at receivers using equation 21 in paper:
        do p = Rx_st,Rx_end
            DBZ_Rx(iter_n,p-Rx_st) = ((EX(i0-1,p,k0,1)-EX(i0-1,p-1,k0,1))/dy(p-1) &
            -(EY(i0,p-1,k0,1)-EY(i0-1,p-1,k0,1))/dx(i0-1))
        end do
        ! The electromagnetic field are used for the next iteration:
        EX(:,:,:,0) = EX(:,:,:,1)
        EY(:,:,:,0) = EY(:,:,:,1)
        EZ(:,:,:,0) = EZ(:,:,:,1)
        HX(:,:,:,0) = HX(:,:,:,1)
        HY(:,:,:,0) = HY(:,:,:,1)
        HZ(:,:,:,0) = HZ(:,:,:,1)
        dt(0) = dt(1)
        iter_n = iter_n+1         
    end do
    print*,"Computation finished."
    !-----------------------------------------------------------------------------
    ! [4.0] Output results (time and dBz/dt at receivers):
    !-----------------------------------------------------------------------------
    deallocate(dx,dy,dz,model_EC,EC_x,EC_y,EC_z,EX,EY,EZ,HX,HY,HZ) 
    ! Output time:
    open(11,file = './Data/Result_time.txt',ACCESS = 'sequential')
    do i = 0,iter_n_max-1
        write(11,*) t_iteration_H(i)
    end do
    close(11)
    ! Output dBz/dt:
    write(n_Rx,*) Rx_end-Rx_st+1
    open(12,file = './Data/Result_dBz.txt',access = 'sequential')
    do i = 0,iter_n_max-1
        write(12,'('//n_Rx//'e20.10e3)') DBZ_Rx(i,0:Rx_end-Rx_st)
    end do
    close(12)
    deallocate(t_iteration_E,t_iteration_H,DBZ_Rx)
    EndSubroutine Sub_Iteration