    !+ Read model and configuration parameters from the input file.
    !
    Subroutine Sub_ReadData
    !
    ! In: inputfile
    ! Out: L_loop, XI, YJ, ZK, dx, dy, dz, i0, j0, k0, Rx_st, Rx_end,
    !      iter_n_max, n_subloop, Alpha, model_EC
    ! Description:
    ! Read model and configuration parameters from the input file.
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
    ! Declarations:
    ! Modules used:
    use Mod_Parameters 
    ! Imported Parameters/Variables with intent:
    ! inputfile       ! Name of the input file.
    ! L_loop          ! Length of Tx loop. Current is 1 A by default.
    ! XI, YJ, ZK      ! Number of cells in the x-, y-, and z-directions.
    ! dx, dy, dz      ! Grid size in the x-, y-, and z-directions.
    ! i0, j0, k0      ! Tx loop center is at the bottom face center
                      ! of the cell (i0,j0,k0).
    ! Rx_st, Rx_end   ! Receivers are at the bottom face center of cells
                      ! from (i0,Rx_st,k0) to (i0,Rx_end,k0).
    ! iter_n_max      ! Maximum number of iterations.
    ! n_subloop       ! Number of subloops.
    ! Alpha           ! Coefficient in the time step size calculation equation.
    ! model_EC        ! Model conductivity.
    implicit none
    ! Subroutine arguments:
    integer::i             ! Temporary loop variables.
    integer::n_submodel    ! Number of model subareas.
    real(8)::submodel(0:6) ! Temporary array used to read model parameters.
    !- End of header --------------------------------------------------------------

    open(111,file = inputfile,status = 'old') ! Open the input file.
    read(111,*)L_loop      ! Read length of Tx loop (m). Current is 1 A by default.
    read(111,*)XI,YJ,ZK    ! Read number of cells in the x-, y-, and z-directions.
    allocate (dx(0:XI-1))
    allocate (dy(0:YJ-1))
    allocate (dz(0:ZK-1))
    read(111,*) dx(0:XI-1) ! Read grid size in the x direction (m).
    read(111,*) dy(0:YJ-1) ! Read grid size in the y direction (m).
    read(111,*) dz(0:ZK-1) ! Read grid size in the z direction (m).
    read(111,*) i0,j0,k0   ! Read position of Tx loop. Tx loop center is at 
                           ! the bottom face center of the cell (i0,j0,k0).
    read(111,*)Rx_st,Rx_end! Read positions of the start and the end receivers. 
                           ! Receivers are at the bottom face center of cells 
                           ! from (i0,Rx_st,k0) to (i0,Rx_end,k0).
    read(111,*)iter_n_max  ! Read maximum number of iterations.
    read(111,*)n_subloop   ! Read number of subloops.
    read(111,*) Alpha      ! Read coefficient in the time step size calculation equation.
    Call Sub_AllocateParameters ! Allocate parameters & variables used in this software.
    read(111,*) n_submodel ! Read number of model subareas.
    do i = 1,n_submodel    ! Read ranges and conductivities of subareas (S/m).
        read(111,*) submodel(0:6)
        model_EC(submodel(0)-1:submodel(1)-1,submodel(2)-1:submodel(3)-1,  &
        submodel(4)-1:submodel(5)-1) = submodel(6)
    end do      
    close(111)
    EndSubroutine Sub_ReadData