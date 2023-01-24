    !+ Global parameters and variables used in this software.
    !
    Module Mod_Parameters
    !
    ! Description:
    ! Global parameters and variables used in this software.
    !
    ! Current Code Owner: <Fei Li and Jiulong Cheng>
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
    ! Global Parameters:
    parameter(pi = 3.141592653589793d0)
    parameter(Permeability_vac = 1.256637061435917d-6) ! Magnetic permeability of vacuum.
    
    ! Global Scalars:
    character(len = 256)::inputfile          ! Name of the input file.
    integer::XI,YJ,ZK       ! Number of cells in the x-, y-, and z-directions.
    integer::i0,j0,k0       ! Tx loop center is at the bottom face center
                            ! of the cell (i0,j0,k0).
    integer::Rx_st,Rx_end   ! Receivers are at the bottom face center of cells
                            ! from (i0,Rx_st,k0) to (i0,Rx_end,k0).
    integer::iter_n_max                      ! Maximum number of iterations.
    real(8)::L_loop                          ! Length of Tx loop.
    real(8)::n_subloop                       ! Number of subloops.
    real(8)::L_subloop                       ! Length of subloops.
    real(8)::subloop_x,subloop_y             ! Coordinate of the subloop center.
    real(8)::t0_E                            ! Initial time of electric field.
    real(8)::t0_H                            ! Initial time of magnetic field.
    real(8)::Alpha      ! Coefficient in the time step size calculation equation.
    
    ! Global Arrays:
    real(8),allocatable::dx(:)               ! Grid size in the x direction.
    real(8),allocatable::dy(:)               ! Grid size in the y direction.
    real(8),allocatable::dz(:)               ! Grid size in the z direction.
    real(8),allocatable::model_EC(:,:,:)     ! Model conductivity.
    real(8),allocatable::EC_x(:,:,:)         ! Conductivity of nodes.
    real(8),allocatable::EC_y(:,:,:)         ! ditto.
    real(8),allocatable::EC_z(:,:,:)         ! ditto.
    real(8),allocatable::t_iteration_E(:)    ! Time of electric field.
    real(8),allocatable::t_iteration_H(:)    ! Time of magnetic field.
    real(8)::dt(0:1)                         ! Time step size.
    real(8),allocatable::EX(:,:,:,:)         ! Ex.
    real(8),allocatable::EY(:,:,:,:)         ! Ey.
    real(8),allocatable::EZ(:,:,:,:)         ! Ez.
    real(8),allocatable::HX(:,:,:,:)         ! Hx.
    real(8),allocatable::HY(:,:,:,:)         ! Hy.
    real(8),allocatable::HZ(:,:,:,:)         ! Hz.
    real(8),allocatable::EX_subloop(:,:,:,:) ! Subloop initial field Ex.
    real(8),allocatable::EY_subloop(:,:,:,:) ! Subloop initial field Ey.
    real(8),allocatable::EZ_subloop(:,:,:,:) ! Subloop initial field Ez.
    real(8),allocatable::HX_subloop(:,:,:,:) ! Subloop initial field Hx.
    real(8),allocatable::HY_subloop(:,:,:,:) ! Subloop initial field Hy.
    real(8),allocatable::HZ_subloop(:,:,:,:) ! Subloop initial field Hz.
    real(8),allocatable::DBZ_Rx(:,:)         ! dBz/dt at receivers.

    EndModule Mod_Parameters

    !- End of module header