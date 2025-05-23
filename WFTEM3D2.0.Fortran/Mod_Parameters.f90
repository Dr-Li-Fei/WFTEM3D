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
    ! 1.0      01/10/21  Original code. Fei Li
    !- End of module header --------------------------------------------------------------

    ! Global Parameters:
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: mu_0 = 1.256637061435917d-6 ! Permeability of vacuum.
    
    ! Global Scalars:
    character(len = 256)::inputfile          ! Name of the input file.
    integer::XI,YJ,ZK       ! Number of cells in the x-, y-, and z-directions.                      
    integer::iter_n_pri                      ! Iteration number of primary field.
    integer::iter_n_sec                      ! Iteration number of secondary field.
    real(8)::Alpha_pri      ! Coefficient of time step size for primary field.
    real(8)::Alpha_sec      ! Coefficient of time step size for secondary field.

    ! Global Arrays:
    real(8),allocatable::dx(:)               ! Grid size in the x direction.
    real(8),allocatable::dy(:)               ! Grid size in the y direction.
    real(8),allocatable::dz(:)               ! Grid size in the z direction.
    integer::Tx(0:5)                         ! Position of Tx.
    integer::Rx(0:8)                         ! Position of Rx.
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
    real(8),allocatable::DBZ_Rx(:,:)         ! dBz/dt at receivers.

    EndModule Mod_Parameters