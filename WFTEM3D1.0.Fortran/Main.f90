    !+ A quick and flexible 3D transient electromagnetic modeling software.
    !
    Program WFTEM3D

    ! Description:
    ! WFTEM3D is a quick, flexible and extensible 3D finite-difference 
    ! transient electromagnetic modeling open-source software. It is 
    ! applicable to the tunnel model, the non-flat topography model, 
    ! the multi-scale model, etc.
    !
    ! Method:
    ! First the scheme calculates the approximate initial field excited by
    ! whole-space magnetic dipole sources at an initial time after the current 
    ! is switched off. Then the scheme steps Maxwell¡¯s equations in time 
    ! using a staggered grid and a modified DuFort-Frankel method. See the 
    ! paper "3D Finite-difference Transient Electromagnetic Modeling with 
    ! Whole-space Initial Field" for detail.
    !
    ! Input files:
    ! ".dat" file,
    ! e.g., "Example1_Conductive brick in a half-space.dat"
    !     , "Example2_Complex conductor at a vertical contact.dat".
    ! In which routine they are read: Sub_ReadData.
    !
    ! Output files:
    ! "Result_time.txt", "Result_dBz.txt" and "Run_time.txt".
    ! In which routine they are written: Sub_Iteration.
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
    ! character(len = 256)::inputfile        ! Name of the input file.
    implicit none
    ! Subroutine arguments:
    real(8)::CPU_time_begin,CPU_time_end     ! The begin and end of run-time.
    !- End of header --------------------------------------------------------------

    call CPU_time(CPU_time_begin)
    !-----------------------------------------------------------------------------
    ! [1.0] Read model and configuration parameters from the input file:
    !-----------------------------------------------------------------------------    
    inputfile = './Data/Example1_Conductive brick in a half-space.dat'
    !inputfile = '../Data/Example2_Complex conductor at a vertical contact.dat'
    call Sub_ReadData
    !-----------------------------------------------------------------------------
    ! [2.0] Calculate initial field:
    !-----------------------------------------------------------------------------    
    call Sub_InitialField
    !-----------------------------------------------------------------------------
    ! [3.0] Iteration and output results:
    !-----------------------------------------------------------------------------    
    call Sub_Iteration
    call CPU_time(CPU_time_end)
    ! Output run-time:
    open(13,file = './Data/Run_time.txt',access = 'sequential')
    write(13,*) 'Computation finished. Run-time is ',CPU_time_end-CPU_time_begin,'s.'
    close(13)
    End
    