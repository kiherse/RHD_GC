program Winberger
    implicit none
    real(8), dimension(:), allocatable :: ic_radius, ic_density, ic_pressure, ic_spec_int_energy
    real(8), allocatable :: rad, ic_density_interp, ic_spec_int_energy_interp, ic_pressure_interp
    integer :: n, i
    real(8) :: rb

    !---------------------------------------------------- WINBERGER

    ! Profile initial conditions
    call read_data(ic_radius, ic_density, ic_pressure, ic_spec_int_energy)

    ! Convert radius to kpc
    ic_radius = ic_radius / 3.085678d+21

    ! RATPENAT RESOLUTION
    !rad = linspace(1.0d0, 256.0d0, 512)
    rb = 222 !pc
    rad = 10

    ! INTERPOLATION

    allocate(ic_density_interp, ic_spec_int_energy_interp, ic_pressure_interp)

    ic_density_interp = interpolate(rad * rb / 1.0d3, ic_radius, ic_density)
    ic_spec_int_energy_interp = interpolate(rad * rb / 1.0d3, ic_radius, ic_spec_int_energy)
    ic_pressure_interp = interpolate(rad * rb / 1.0d3, ic_radius, ic_pressure)

    print*,'Original radius', ic_radius(1:5)
    print*,'Original density', ic_density(1:5)
    print*,'Interpolation radius', rad * rb / 1.0d3
    print*,'Interpolated density', ic_density_interp

contains 

    subroutine read_data(ic_radius, ic_density, ic_pressure, ic_spec_int_energy)
        real(8), dimension(:), allocatable :: ic_radius, ic_density, ic_pressure, ic_spec_int_energy
        integer :: unit_number, n

        open(unit=unit_number, file='../perseus_ics.txt', status='unknown', action='read')
        n = count_lines(unit_number) - 1
        allocate(ic_radius(n), ic_density(n), ic_pressure(n), ic_spec_int_energy(n))
        read(unit_number,*) ! Skip header (first line)
        do i = 1, n
            read(unit_number, *) ic_radius(i), ic_density(i), ic_pressure(i), ic_spec_int_energy(i)        
        end do

        close(unit_number)
    end subroutine read_data

    function count_lines(unit_number) result(num_lines)
        integer, intent(in) :: unit_number
        integer :: num_lines, iostat
        character(len=100) :: line
    
        num_lines = 0
    
        do
            read(unit_number, '(a)', iostat=iostat) line
            if (iostat /= 0) exit
            num_lines = num_lines + 1
        end do
    
        rewind(unit_number)
    end function count_lines
    
    function interpolate(x, xp, fp) result(fp_interp)
        real(8), intent(in) :: x
        real(8), dimension(:), intent(in) :: xp, fp
        real(8) :: fp_interp
        integer :: j
    
        j = 1
        do while (xp(j) < x .and. j < size(xp))
            j = j + 1
        end do
        if (j == 1) then
            fp_interp = fp(1)
        else if (j == size(xp)) then
            fp_interp = fp(size(xp))
        else
            fp_interp = fp(j-1) + (x - xp(j-1)) * (fp(j) - fp(j-1)) / (xp(j) - xp(j-1))
        end if

    end function interpolate

end program Winberger
