    !+ Iteration and output results.
    !
    Subroutine Sub_PrimaryField
    !
    !-----------------------------------------------------------------------------
    ! In: XI,YJ,ZK,dx,dy,dz,Tx,iter_n_pri,alpha_pri
    ! Out: EX,EY,EZ,HX,HY,HZ
    ! Description:
    ! Calculatie primary field.
    !
    ! Method:
    ! Wire is modeled as volume current.First calculate the primary field using 
    ! a whole-space homogeneous model, then calculate the secondary field using  
    ! the true model. (Relevant paper will be provided in subsequent updates.)
    !
    ! Current Code Owner: <Fei Li and Jiulong Cheng>
    !
    ! History:
    ! Version    Date    Comment
    ! -------    ----    -------
    ! 1.0      01/10/21  Original code. <Fei Li and Jiulong Cheng>
    ! 2.0      01/03/25  Magnetic dipole-based -> volume current-based. Fei Li
    !
    ! Declarations:
    ! Modules used:
    use Mod_Parameters 
    ! Imported Parameters/Variables with intent:
    ! XI, YJ, ZK           ! Number of cells in the x-, y-, and z-directions.
    ! dx, dy, dz           ! Grid size in the x-, y-, and z-directions.
    ! Tx                   ! Tx loop diagonal points / ground wire endpoints.
    ! iter_n_pri           ! Iteration number of primary field.
    ! alpha_pri            ! Coefficient of time step size for primary field
    ! EX(:,:,:,:)          ! Ex.
    ! EY(:,:,:,:)          ! Ey.
    ! EZ(:,:,:,:)          ! Ez.
    ! HX(:,:,:,:)          ! Hx.
    ! HY(:,:,:,:)          ! Hy.
    ! HZ(:,:,:,:)          ! Hz.
    ! mu_0                 ! Magnetic permeability of vacuum.
    ! t_iteration_E(:)     ! Time of electric field.
    ! t_iteration_H(:)     ! Time of magnetic field.
    ! EC                   ! Conductivity of homogeneous model.   
    ! dt(0:1)              ! Time step size.
    implicit none
    ! Subroutine arguments:
    integer::i,j,k,iter_n                 ! Temporary loop control variables.
    integer::Tx_x1,Tx_y1,Tx_x2,Tx_y2,Tx_z ! Position of Tx.
    integer::EC                           ! Conductivity of homogeneous model.
    integer::current                      ! Current.
    real(8)::gamma                        ! Artificial dielectric constant.
    !- End of header --------------------------------------------------------------

    print*,"Calculating primary field..."   
    !-----------------------------------------------------------------------------
    ! [1.0] Model and configuration initialization:
    !----------------------------------------------------------------------------- 
    allocate (t_iteration_E(0:iter_n_pri-1))
    allocate (t_iteration_H(0:iter_n_pri-1))
    allocate (EX(0:XI-1,0:YJ,0:ZK,0:1))
    allocate (EY(0:XI,0:YJ-1,0:ZK,0:1))
    allocate (EZ(0:XI,0:YJ,0:ZK-1,0:1))
    allocate (HX(0:XI,0:YJ-1,0:ZK-1,0:1))
    allocate (HY(0:XI-1,0:YJ,0:ZK-1,0:1))
    allocate (HZ(0:XI-1,0:YJ-1,0:ZK,0:1))
    EX = 0.0d0; EY = 0.0d0; EZ = 0.0d0; HX = 0.0d0; HY = 0.0d0; HZ = 0.0d0
    EC = 10                            ! Conductivity of homogeneous model.
    current=1                          ! Current is 1 A by default
    Tx_x1=Tx(0);Tx_y1=Tx(1);Tx_z=Tx(2) ! Position of Tx.
    Tx_x2=Tx(3);Tx_y2=Tx(4)            ! Position of Tx.
    t_iteration_E(0) = 0               ! Initial time of electric field.
    ! Calculate the time step size using equation 28 in paper:
    dt(0) =2.26*mu_0*EC*dz(Tx_z-1)**2
    ! Initial time of magnetic field:
    t_iteration_H(0) = t_iteration_E(0)+dt(0)/2
    !-----------------------------------------------------------------------------
    ! [2.0] Calculate primary field:
    !-----------------------------------------------------------------------------
    iter_n = 1
    do while (iter_n <= iter_n_pri-1)
        !print*, iter_n
        ! Calculate the electric field time of the iter_nth iteration:
        t_iteration_E(iter_n) = t_iteration_E(iter_n-1)+dt(0)
        ! Calculate the time step size using equation 23 in paper:       
        dt(1) = alpha_pri*dz(Tx_z-1)*sqrt(mu_0* &
        EC*t_iteration_E(iter_n)/6.0D0)
        ! Calculate the magnetic field time of the iter_nth iteration:
        t_iteration_H(iter_n) = t_iteration_H(iter_n-1)+(dt(0)+dt(1))/2.0d0
        ! Calculate artificial dielectric constant using equation 22 in paper:
        gamma = (4.0d0/mu_0)*(dt(0)/dz(Tx_z-1))**2
        ! Update the value of Ex using equation 12 in paper:
        do k = 1,ZK-1
            do j = 1,YJ-1 
                do i = 0,XI-1
                    EX(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EX(i,j,k,0)  &                          
                    +4.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  & 
                    *((HZ(i,j,k,0)-HZ(i,j-1,k,0))/(dy(j-1)+dy(j)) &                          
                    -(HY(i,j,k,0)-HY(i,j,k-1,0))/(dz(k-1)+dz(k))) 
                end do
            end do   
        end do
        ! Wire along x-direction if Tx_x1 .ne. Tx_x2:
        if (Tx_x1 .ne. Tx_x2) then
            DO i = Tx_x1-1, Tx_x2-1
                j = Tx_y1-1
                k = Tx_z
                    EX(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EX(i,j,k,0)  &                          
                    +2.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  & 
                    *((HZ(i,j,k,0)-HZ(i,j-1,k,0))/((dy(j-1)+dy(j))/2.0d0) &                          
                    -(HY(i,j,k,0)-HY(i,j,k-1,0))/((dz(k-1)+dz(k))/2.0d0)  &
                    -current/((dy(j-1)+dy(j))*(dz(k-1)+dz(k))/4.0D0)) 
            end do
        end if
        ! Anothor wire along x-direction if Tx_x1 .ne. Tx_x2 and Tx_y1 .ne. Tx_y2:
        if (Tx_x1 .ne. Tx_x2 .and. Tx_y1 .ne. Tx_y2) then
            do i = Tx_x1-1, Tx_x2-1
                j = Tx_y2
                k = Tx_z
                    EX(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EX(i,j,k,0)  &                          
                    +2.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  & 
                    *((HZ(i,j,k,0)-HZ(i,j-1,k,0))/((dy(j-1)+dy(j))/2.0d0) &                          
                    -(HY(i,j,k,0)-HY(i,j,k-1,0))/((dz(k-1)+dz(k))/2.0d0)  &
                    -current*(-1.0d0)/((dy(j-1)+dy(j))*(dz(k-1)+dz(k))/4.0D0)) 
            end do
        end if
        ! Update the value of Ey using equation 14 in paper:
        do k = 1,ZK-1
            do j = 0,YJ-1 
                do i = 1,XI-1             
                    EY(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EY(i,j,k,0)  &
                    +4.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  &
                    *((HX(i,j,k,0)-HX(i,j,k-1,0))/(dz(k-1)+dz(k)) &            
                    -(HZ(i,j,k,0)-HZ(i-1,j,k,0))/(dx(i-1)+dx(i)))         
                end do
            end do   
        end do
        ! Wire along y-direction if Tx_y1 .ne. Tx_y2:
        if (Tx_y1 .ne. Tx_y2) then
            do j = Tx_y1-1, Tx_y2-1
                i = Tx_x1-1
                k = Tx_z
                    EY(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EY(i,j,k,0)  &
                    +2.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  &
                    *((HX(i,j,k,0)-HX(i,j,k-1,0))/((dz(k-1)+dz(k))/2.0d0) &            
                    -(HZ(i,j,k,0)-HZ(i-1,j,k,0))/((dx(i-1)+dx(i))/2.0d0)  &
                    -current*(-1.0d0)/((dz(k-1)+dz(k))*(dx(i-1)+dx(i))/4.0d0))
            end do
        end if
        ! Anothor wire along y-direction if Tx_x1 .ne. Tx_x2 and Tx_y1 .ne. Tx_y2:
        if (Tx_x1 .ne. Tx_x2 .and. Tx_y1 .ne. Tx_y2) then
            do j = Tx_y1-1, Tx_y2-1
                i = Tx_x2
                k = Tx_z
                    EY(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EY(i,j,k,0)  &
                    +2.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  &
                    *((HX(i,j,k,0)-HX(i,j,k-1,0))/((dz(k-1)+dz(k))/2.0d0) &            
                    -(HZ(i,j,k,0)-HZ(i-1,j,k,0))/((dx(i-1)+dx(i))/2.0d0)  &
                    -current/((dz(k-1)+dz(k))*(dx(i-1)+dx(i))/4.0d0))
            end do
        end if
        ! Update the value of Ez using equation 16 in paper:
        do k = 0,ZK-1     
            do j = 1,YJ-1      
                do i = 1,XI-1             
                    EZ(i,j,k,1) = (2.0d0*gamma-EC*dt(0)) &
                    /(2.0d0*gamma+EC*dt(0))*EZ(i,j,k,0)  &
                    +4.0d0*dt(0)/(2.0d0*gamma+EC*dt(0))  &
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
            print '(I6, A, F5.1, A)', iter_n, " (", 100.0*iter_n/iter_n_pri, "%)"
        end if         
    end do
    deallocate(t_iteration_E,t_iteration_H)  
    print*,"Primary field calculation finished."

    EndSubroutine Sub_PrimaryField