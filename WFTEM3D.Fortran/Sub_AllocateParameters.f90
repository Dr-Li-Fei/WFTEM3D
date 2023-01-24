    !+ Allocate parameters and variables used in this software.
    !
    Subroutine Sub_AllocateParameters
    !
    ! Description:
    ! Allocate parameters and variables used in this software.
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
    ! model_EC(:,:,:)       ! Model conductivity.
    ! EC_x(:,:,:)           ! Conductivity of nodes.
    ! EC_y(:,:,:)           ! ditto.
    ! EC_z(:,:,:)           ! ditto .  
    ! t_iteration_E(:)      ! Time of electric field.
    ! t_iteration_H(:)      ! Time of magnetic field.
    ! EX(:,:,:,:)           ! Ex.
    ! EY(:,:,:,:)           ! Ey.
    ! EZ(:,:,:,:)           ! Ez.
    ! HX(:,:,:,:)           ! Hx.
    ! HY(:,:,:,:)           ! Hy.
    ! HZ(:,:,:,:)           ! Hz.
    ! EX_subloop(:,:,:,:)   ! Subloop initial field Ex.
    ! EY_subloop(:,:,:,:)   ! Subloop initial field Ey.
    ! EZ_subloop(:,:,:,:)   ! Subloop initial field Ez.
    ! HX_subloop(:,:,:,:)   ! Subloop initial field Hx.
    ! HY_subloop(:,:,:,:)   ! Subloop initial field Hy.
    ! HZ_subloop(:,:,:,:)   ! Subloop initial field Hz.
    ! DBZ_Rx(:,:)           ! dBz/dt at receivers. 
    implicit none
    !- End of header --------------------------------------------------------------
    
    allocate (model_EC(0:XI-1,0:YJ-1,0:ZK-1))  
    allocate (EC_x(0:XI-1,0:YJ,0:ZK))
    allocate (EC_y(0:XI,0:YJ-1,0:ZK))
    allocate (EC_z(0:XI,0:YJ,0:ZK-1))
    allocate (t_iteration_E(0:iter_n_max-1))
    allocate (t_iteration_H(0:iter_n_max-1))
    allocate (EX(0:XI-1,0:YJ,0:ZK,0:1))
    allocate (EY(0:XI,0:YJ-1,0:ZK,0:1))
    allocate (EZ(0:XI,0:YJ,0:ZK-1,0:1))
    allocate (HX(0:XI,0:YJ-1,0:ZK-1,0:1))
    allocate (HY(0:XI-1,0:YJ,0:ZK-1,0:1))
    allocate (HZ(0:XI-1,0:YJ-1,0:ZK,0:1))
    allocate (EX_subloop(0:XI-1,0:YJ,0:ZK,0:1))
    allocate (EY_subloop(0:XI,0:YJ-1,0:ZK,0:1))
    allocate (EZ_subloop(0:XI,0:YJ,0:ZK-1,0:1))
    allocate (HX_subloop(0:XI,0:YJ-1,0:ZK-1,0:1))
    allocate (HY_subloop(0:XI-1,0:YJ,0:ZK-1,0:1))
    allocate (HZ_subloop(0:XI-1,0:YJ-1,0:ZK,0:1))
    allocate (DBZ_Rx(0:iter_n_max-1,0:Rx_end-Rx_st))
    EX = 0.0d0
    EY = 0.0d0
    EZ = 0.0d0
    HX = 0.0d0
    HY = 0.0d0
    HZ = 0.0d0
    EndSubroutine Sub_AllocateParameters
    