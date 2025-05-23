    !+ A quick and flexible 3D transient electromagnetic modeling software.
    !
    Program WFTEM3D

    ! Description:
    ! A quick, flexible, and extensible 3D TEM modeling open-source software. 
    ! Supports wire/loop sources, half/whole-space, and ground/airborne/marine/ 
    ! tunnel/borehole scenarios. 
    !
    ! Method:
    ! Wire is modeled as volume current.First calculate the primary field using
    ! a whole-space homogeneous model, then calculate the secondary field using 
    ! the true model. (Relevant paper will be provided in subsequent updates.)
    !
    ! Input files:
    ! ".txt" file or ".dat" file,
    ! e.g., "Example_conductive_brick_in_a_half-space.txt"
    ! In which routine they are read: Sub_ReadData.
    !
    ! Output files:
    ! "Result_time.txt" and "Result_dBz.txt".
    ! In which routine they are written: Sub_SecondaryField.
    !
    ! Current Code Owner: <Fei Li and Jiulong Cheng>
    !
    ! History:
    ! Version    Date    Comment
    ! -------    ----    -------
    ! 1.0      01/10/21  Original code. Fei Li
    ! 2.0      01/03/25  Add wire source support(volume current-based). Fei Li
    !
    ! Declarations:
    ! Modules used:
    use Mod_Parameters
    ! Imported Parameters/Variables with intent:
    ! character(len = 256)::inputfile        ! Name of the input file.
    ! Subroutine arguments:
    real(8)::CPU_time_begin,CPU_time_end     ! The begin and end of run-time.
    !- End of header --------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [1.0] Read model and configuration parameters from the input file:
    !-----------------------------------------------------------------------------    
    write(*,*) "Enter input filename (e.g., input.txt or input.dat):"
    read(*,*) inputfile
    call Sub_ReadData
    call CPU_time(CPU_time_begin)
    !-----------------------------------------------------------------------------
    ! [2.0] Calculate primary field:
    !----------------------------------------------------------------------------- 
    call Sub_PrimaryField
    !-----------------------------------------------------------------------------
    ! [3.0] Calculate secondary field:
    !----------------------------------------------------------------------------- 
    call Sub_SecondaryField
    call CPU_time(CPU_time_end)
    ! Output run-time:
    open(13,file = 'Run_time.txt',access = 'sequential')
    write(13,*) 'Computation finished. Run-time is ',CPU_time_end-CPU_time_begin,'s.'
    close(13)
    read(*,*)
    End