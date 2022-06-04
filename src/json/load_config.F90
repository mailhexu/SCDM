    program example1

        use json_module

        implicit none

        type(json_file) :: json
        logical :: found
        integer :: int_var, missing_int_var
        real(8) :: float_var
        integer, allocatable :: int_list_var(:)
        real(8), allocatable :: float_list_var(:)
        CHARACTER(len=:), allocatable :: string_var

        ! initialize the class
        call json%initialize()

        ! read the file
        call json%load(filename = 'input.json')

        ! print the file to the console
        call json%print()

        ! extract data from the file
        ! [found can be used to check if the data was really there]
        call json%get('int_var', int_var, found)
        if ( .not. found ) stop 1
        call json%get('float_var', float_var, found)
        if ( .not. found ) stop 1

        call json%get('int_list_var', int_list_var, found)
        if ( .not. found ) stop 1

        call json%get('float_list_var', float_list_var, found)
        if ( .not. found ) stop 1

        call json%get('string_var', string_var, found)
        if ( .not. found ) stop 1


        ! missing value: use default
        
        call json%get('missing_int_var', missing_int_var, found)
        if (.not. found) missing_int_var=-99

        !call json%get('data(1).number', k, found)
        !if ( .not. found ) stop 1

        print *, "int_var: ", int_var
        print *, "missing_int_var: ", int_var
        print *, "float_var: ", float_var
        print *, "int_list_var: ", int_list_var
        print *, "float_list_var: ", float_list_var
        print *, "string_var: ", string_var

        ! clean up
        
        deallocate(int_list_var)
        deallocate(float_list_var)
        deallocate(string_var)

        call json%destroy()
        if (json%failed()) stop 1

    end program example1
