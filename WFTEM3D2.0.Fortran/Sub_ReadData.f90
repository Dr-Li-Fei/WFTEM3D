    !+ Read model and configuration parameters from the input file.
    !
    Subroutine Sub_ReadData
    !
    ! In: inputfile
    ! Out: XI,YJ,ZK,dx,dy,dz,Tx,Rx,iter_n_pri,iter_n_sec,alpha_pri,
    !      alpha_sec,model_EC
    ! Description:
    ! Read model and configuration parameters from the input file.
    !
    ! Current Code Owner: <Fei Li and Jiulong Cheng>
    !
    ! History:
    ! Version    Date    Comment
    ! -------    ----    -------
    ! 1.0      01/10/21  Original code. Fei Li
    ! 2.0      01/03/25  Add support for annotated input file. Fei Li
    !
    ! Declarations:
    ! Modules used:
    use Mod_Parameters
    ! Imported Parameters/Variables with intent:
    ! inputfile      ! Name of the input file.
    ! XI, YJ, ZK     ! Number of cells in the x-, y-, and z-directions.
    ! dx, dy, dz     ! Grid size in the x-, y-, and z-directions.
    ! Tx             ! Tx loop diagonal points / ground wire endpoints.
    ! Rx             ! Rx points: start/end/interval.
    ! iter_n_pri     ! Iteration number of primary field.
    ! iter_n_sec     ! Iteration number of secondary field.
    ! alpha_pri      ! Coefficient of time step size for primary field.
    ! alpha_sec      ! Coefficient of time step size for secondary field.
    ! model_EC       ! Model conductivity.
    implicit none
    ! Subroutine arguments:
    character(len=25600)::line  ! Temporary variable for storing each line.
    logical::file_exists    ! Check file existence before opening.
    integer::i              ! Temporary loop variables.
    integer::n_submodel     ! Number of subdomains in the model.
    integer::submodel(0:5)  ! Index ranges for each subdomain.
    real(8)::submodel_EC    ! Electrical conductivity for each subdomain.
    !- End of header --------------------------------------------------------------
    
    inquire(file=inputfile, exist=file_exists)! Check file existence.
    if (.not. file_exists) then
        write(*,*) "ERROR: Input file '"//trim(inputfile)//"' not found!"
        read(*,*)
    end if

    open(111,file = inputfile,status = 'old') ! Open the input file.
    call skip_comments(111, line)
    read(line,*) XI,YJ,ZK   ! Read number of cells in the x-, y-, and z-directions.
    allocate (dx(0:XI-1), dy(0:YJ-1), dz(0:ZK-1))
    call skip_comments(111, line)
    read(line,*) dx(0:XI-1) ! Read grid size in the x direction (m).
    call skip_comments(111, line)
    read(line,*) dy(0:YJ-1) ! Read grid size in the y direction (m).
    call skip_comments(111, line)
    read(line,*) dz(0:ZK-1) ! Read grid size in the z direction (m).
    call skip_comments(111, line)
    read(line,*) Tx(0:2)    ! Read position of Tx.
    call skip_comments(111, line)
    read(line,*) Tx(3:5)    ! Read position of Tx.
    call skip_comments(111, line)
    read(line,*) Rx(0:2)    ! Read position of Rx.
    call skip_comments(111, line)
    read(line,*) Rx(3:5)    ! Read position of Rx.
    call skip_comments(111, line)
    read(line,*) Rx(6:8)    ! Read position of Rx.
    call skip_comments(111, line)
    read(line,*) iter_n_pri,iter_n_sec ! Read iteration number of primary field.
    call skip_comments(111, line)
    read(line,*) Alpha_pri,Alpha_sec   ! Read iteration number of secondary field.
    allocate (model_EC(0:XI-1,0:YJ-1,0:ZK-1))
    call skip_comments(111, line)
    read(line,*) n_submodel ! Read number of model subdomains.
    do i = 1,n_submodel     ! Read ranges and conductivities of subdomains (S/m).
        call skip_comments(111, line)
        read(line,*) submodel(0:5), submodel_EC
        model_EC(submodel(0)-1:submodel(1)-1,submodel(2)-1:submodel(3)-1,  &
        submodel(4)-1:submodel(5)-1) = submodel_EC
    end do      
    close(111)

    contains ! Subroutine: Skip comment lines (starting with #)
    subroutine skip_comments(file_unit, line)
        integer, intent(in)::file_unit
        character(len=*), intent(out)::line
        do
            read(file_unit, '(a)') line
            line = adjustl(line)
            if (line(1:1) /= '#' .and. line /= '') exit
        end do
    end subroutine skip_comments

    EndSubroutine Sub_ReadData