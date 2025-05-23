    !+ Iteration and output results.
    !
    Subroutine Sub_SecondaryField
    !
    !-----------------------------------------------------------------------------
    ! In: XI,YJ,ZK,dx,dy,dz,Rx,Tx,iter_n_sec,alpha_sec,model_EC,
    !     EX,EY,EZ,HX,HY,HZ
    ! Out: t_iteration_H,DBZ_Rx
    ! Description:
    ! Calculatie secondary field and output results.
    !
    ! Method:
    ! Step Maxwell’s equations in time using a staggered grid and 
    ! a modified DuFort-Frankel method. See 3D Finite-difference Transient 
    ! Electromagnetic Modeling with a Whole-space Initial Field for detail.
    !
    ! Current Code Owner: <Fei Li and Jiulong Cheng>
    !
    ! History:
    ! Version    Date    Comment
    ! -------    ----    -------
    ! 1.0      01/10/21  Original code. Fei Li
    !
    ! Declarations:
    ! Modules used:
    use Mod_Parameters 
    ! Imported Parameters/Variables with intent:
    ! mu_0                ! Magnetic permeability of vacuum.
    ! XI, YJ, ZK          ! Number of cells in the x-, y-, and z-directions.
    ! dx, dy, dz          ! Grid size in the x-, y-, and z-directions.
    ! Tx                  ! Tx loop diagonal points / ground wire endpoints.
    ! Rx                  ! Rx points: start/end/interval.
    ! iter_n_sec          ! Iteration number of secondary field.
    ! alpha_sec           ! Coefficient of time step size for secondary field.
    ! model_EC            ! Model conductivity.
    ! EX(:,:,:,:)         ! Ex.
    ! EY(:,:,:,:)         ! Ey.
    ! EZ(:,:,:,:)         ! Ez.
    ! HX(:,:,:,:)         ! Hx.
    ! HY(:,:,:,:)         ! Hy.
    ! HZ(:,:,:,:)         ! Hz.
    ! DBZ_Rx(:,:)         ! dBz/dt at receivers.
    ! t_iteration_E(:)    ! Time of electric field.
    ! t_iteration_H(:)    ! Time of magnetic field.
    ! EC_x(:,:,:)         ! Conductivity of nodes.
    ! EC_y(:,:,:)         ! ditto.
    ! EC_z(:,:,:)         ! ditto.
    ! dt(0:1)             ! Time step size.
    implicit none
    ! Subroutine arguments:
    integer::i,j,k,iter_n           ! Temporary loop control variables.
    integer::Tx_z                   ! Transmitter depth.
    real(8)::gamma                  ! Artificial dielectric constant.
    real(8)::model_EC_min           ! Minimum conductivity in model.
    integer::Rx_start_x, Rx_end_x, Rx_step_x ! Rx along x-axis:start/end/interval.
    integer::Rx_start_y, Rx_end_y, Rx_step_y ! Rx along y-axis:start/end/interval.
    integer::Rx_start_z, Rx_end_z, Rx_step_z ! Rx along z-axis:start/end/interval.
    integer::num_Rx                          ! Number of receivers.
    integer,allocatable :: R(:,:)            ! Position of each receiver.
    !- End of header --------------------------------------------------------------

    print*,"Calculating secondary field..."
    !-----------------------------------------------------------------------------
    ! [1.0] Model and configuration initialization:
    !-----------------------------------------------------------------------------
    allocate (t_iteration_E(0:iter_n_sec-1))
    allocate (t_iteration_H(0:iter_n_sec-1))   
    Tx_z=Tx(2)                        ! Transmitter depth
    model_EC_min = minval(model_EC)   ! Minimum conductivity in model
    ! Initial time of electric field:
    t_iteration_E(0) = 1.13*mu_0*model_EC_min*dz(Tx_z)**2 
    ! Calculate the time step size using equation 28 in paper:
    dt(0) = Alpha_pri*dz(Tx_z-1)*sqrt(mu_0*model_EC_min*t_iteration_E(0)/6)
    ! Initial time of magnetic field:
    t_iteration_H(0) = t_iteration_E(0)+dt(0)/2
    ! Receiver positions:
    Rx_start_x = Rx(0); Rx_end_x = Rx(1); Rx_step_x = Rx(2)
    Rx_start_y = Rx(3); Rx_end_y = Rx(4); Rx_step_y = Rx(5)
    Rx_start_z = Rx(6); Rx_end_z = Rx(7); Rx_step_z = Rx(8)
    ! Number of receivers:
    num_Rx = ((Rx_end_x - Rx_start_x)/Rx_step_x + 1) * &
         ((Rx_end_y - Rx_start_y)/Rx_step_y + 1) * &
         ((Rx_end_z - Rx_start_z)/Rx_step_z + 1)
    allocate (R(0:num_Rx-1, 0:2))
    allocate (DBZ_Rx(0:iter_n_sec-1, 0:num_Rx-1))
    iter_n = 0
    do k = Rx_start_z-1, Rx_end_z-1, Rx_step_z
        do j = Rx_start_y-1, Rx_end_y-1, Rx_step_y
            do i = Rx_start_x-1, Rx_end_x-1, Rx_step_x
                R(iter_n,:) = [i,j,k] ! Position of each receiver.
                iter_n = iter_n + 1
            end do
        end do
    end do
    ! Calculate dBz/dt at receivers using equation 21 in paper:
    do i = 0, num_Rx-1
        DBZ_Rx(0,i) = ((EX(R(i,0), R(i,1)+1, R(i,2)+1, 0) &
            -EX(R(i,0), R(i,1), R(i,2)+1, 0))/dy(R(i,1))) &
            -((EY(R(i,0)+1, R(i,1), R(i,2)+1, 0)          &
            -EY(R(i,0), R(i,1), R(i,2)+1, 0))/dx(R(i,0)))
    end do
    !-----------------------------------------------------------------------------
    ! [2.0] Calculate conductivity of grid nodes:
    !-----------------------------------------------------------------------------     
    allocate (EC_x(0:XI-1,0:YJ,0:ZK))
    allocate (EC_y(0:XI,0:YJ-1,0:ZK))
    allocate (EC_z(0:XI,0:YJ,0:ZK-1))
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
    ! [3.0] Calculate secondary field:
    !-----------------------------------------------------------------------------
    iter_n = 1
    do while (iter_n <= iter_n_sec-1)
        !print*, iter_n
        ! Calculate the electric field time of the iter_nth iteration:
        t_iteration_E(iter_n) = t_iteration_E(iter_n-1)+dt(0)
        ! Calculate the time step size using equation 23 in paper:       
        dt(1) = Alpha_sec*dz(Tx_z-1)*sqrt(mu_0* &
        model_EC_min*t_iteration_E(iter_n)/6.0D0)
        ! Calculate the magnetic field time of the iter_nth iteration:
        t_iteration_H(iter_n) = t_iteration_H(iter_n-1)+(dt(0)+dt(1))/2.0d0
        ! Calculate artificial dielectric constant using equation 22 in paper:
        gamma = (4.0d0/mu_0)*(dt(0)/dz(Tx_z-1))**2
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
        EX(0:XI-1, 0:YJ, 0, 1) = 0.0d0
        EX(0:XI-1, 0:YJ, ZK, 1) = 0.0d0
        EX(0:XI-1, 0, 0:ZK, 1) = 0.0d0
        EX(0:XI-1, YJ, 0:ZK, 1) = 0.0d0
        EY(0:XI, 0:YJ-1, 0, 1) = 0.0d0
        EY(0:XI, 0:YJ-1, ZK, 1) = 0.0d0
        EY(0, 0:YJ-1, 0:ZK, 1) = 0.0d0
        EY(XI, 0:YJ-1, 0:ZK, 1) = 0.0d0
        EZ(0:XI, 0, 0:ZK-1, 1) = 0.0d0
        EZ(0:XI, YJ, 0:ZK-1, 1) = 0.0d0
        EZ(0, 0:YJ, 0:ZK-1, 1) = 0.0d0
        EZ(XI, 0:YJ, 0:ZK-1, 1) = 0.0d0
        ! Update the value of Hx using equation 18 in paper:
        do k = 0,ZK-1
            do j = 0,YJ-1
                do i = 0,XI
                    HX(i,j,k,1) = HX(i,j,k,0)           &
                    +(dt(0)+dt(1))/(2*mu_0) &
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
                    +(dt(0)+dt(1))/(2*mu_0) &		
                    *((EZ(i+1,j,k,1)-EZ(i,j,k,1))/dx(i) &
                    -(EX(i,j,k+1,1)-EX(i,j,k,1))/dz(k))
                end do
            end do  
        end do
        ! Update the value of Hz using equation 20 in paper, start from 
        ! the bottom boundary and the top boundary of the model where Hz=0:
        HZ(0:XI-1,0:YJ-1,0,1) = 0
        HZ(0:XI-1,0:YJ-1,ZK,1) = 0  
        do k = ZK,Tx_z+1,-1
            do j = 0,YJ-1        
                do i = 0,XI-1
                    HZ(i,j,k-1,1) = HZ(i,j,k,1)                    &
                    +dz(k-1)/dx(i)*(HX(i+1,j,k-1,1)-HX(i,j,k-1,1)) &
                    +dz(k-1)/dy(j)*(HY(i,j+1,k-1,1)-HY(i,j,k-1,1))
                end do
            end do
        end do
        do k = 1,Tx_z-1
            do j = 0,YJ-1        
                do i = 0,XI-1
                    HZ(i,j,k,1) = HZ(i,j,k-1,1)                    &
                    -dz(k-1)/dx(i)*(HX(i+1,j,k-1,1)-HX(i,j,k-1,1)) &
                    -dz(k-1)/dy(j)*(HY(i,j+1,k-1,1)-HY(i,j,k-1,1))
                end do
            end do
        end do
        ! Calculate dBz/dt at receivers using equation 21 in paper:

        do i = 0, num_Rx-1
            DBZ_Rx(iter_n,i) = ((EX(R(i,0), R(i,1)+1, R(i,2)+1, 1) &
                -EX(R(i,0), R(i,1), R(i,2)+1, 1))/dy(R(i,1)))      &
                -((EY(R(i,0)+1, R(i,1), R(i,2)+1, 1)               &
                -EY(R(i,0), R(i,1), R(i,2)+1, 1))/dx(R(i,0)))
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
        ! Print progress message every 500 iterations:
        if (mod(iter_n, 500) == 0) then
            print '(I6, A, F5.1, A)', iter_n, " (", 100.0*iter_n/iter_n_sec, "%)"
        end if          
    end do
    print*,"Secondary field calculation finished."
    !-----------------------------------------------------------------------------
    ! [4.0] Output results (time and dBz/dt at receivers):
    !-----------------------------------------------------------------------------
    deallocate(dx,dy,dz,model_EC,EC_x,EC_y,EC_z,EX,EY,EZ,HX,HY,HZ)
    ! Output time:
    open(11,file = 'Result_time.txt',ACCESS = 'sequential')
    do i = 0,iter_n_sec-1
        write(11,*) t_iteration_H(i)
    end do 
    ! Output dBz/dt:
    open(12, file='Result_dBz.txt')
    do i = 0, iter_n_sec-1
        write(12, '(*(e20.10e3))') DBZ_Rx(i, 0:num_Rx-1)
    end do
    close(12)
    deallocate(t_iteration_E,t_iteration_H,DBZ_Rx)
    EndSubroutine Sub_SecondaryField