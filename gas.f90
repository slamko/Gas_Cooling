
program main
  use, intrinsic :: iso_c_binding, only: c_null_char
  use :: raylib
  implicit none

  external :: sgesv
  integer, parameter :: SCREEN_WIDTH  = 800
  integer, parameter :: SCREEN_HEIGHT = 450
  real, parameter :: MOLEKULA_RADIUS = 6
  real, parameter :: ELASTIC_COEFF = 0.95
  
  type :: dvector
     double precision :: x
     double precision :: y
  end type dvector
 
  type :: particle
     type (vector2_type) :: pos

     type (vector2_type) :: speed
  end type particle
  
  type(particle), dimension (16) :: particles
  integer, parameter :: FPS = 60
  real :: dt
  dt = 1.0 / FPS
 
  print *, 'Hello world'
  
  call random_seed()
  call init_particles(particles, size(particles))
  call print_particles(particles, size(particles))

  call init_window(SCREEN_WIDTH, SCREEN_HEIGHT, 'Fortran raylib' // c_null_char)
  call set_target_fps(fps)
  

  do while (.not. window_should_close())
     call begin_drawing()
     call clear_background(RAYWHITE)
     call draw_particles(particles, size(particles))
     call update_particles(particles, size(particles))
     call end_drawing()

     print *, calc_temp(particles, size(particles))
  end do

  call close_window()
  
contains
  subroutine init_particles (particles, length)
    type(particle), dimension(*) :: particles
    integer :: length
    
    integer :: i
    real :: x (length)
    real :: y (length)
    double precision :: speed (length * 2)
    real :: block_x, block_y

    block_y = SCREEN_HEIGHT / length
    block_x = SCREEN_WIDTH / length

    call random_number(x)
    call random_number(y)
    call random_number(speed)

    do i = 1, length
       particles(i)%pos = vector2_type (block_x * (i - 1) + (block_x / 2), block_y * (i - 1) + (block_y / 2))

       particles(i)%speed = vector2_type ((speed (i * 2) * 2 - 1) * 100, (speed ((i * 2) + 1) * 2 - 1) * 100)
    end do
    
  end subroutine init_particles

  function vsub (v1, v2) result(v3)
    type (vector2_type), intent(in) :: v1, v2
    type (vector2_type) :: v3

    v3%x = v1%x - v2%x
    v3%y = v1%y - v2%y
  end function 

  function vadd (v1, v2) result(v3)
    type (vector2_type), intent(in) :: v1, v2
    type (vector2_type) :: v3

    v3%x = v1%x + v2%x
    v3%y = v1%y + v2%y
  end function 

  function vdot (v1, v2) result(res)
    type (vector2_type), intent(in) :: v1, v2
    real :: res

    res = v1%x * v2%x + v1%y * v2%y
  end function 

  function vscale (v1, fact) result(vec)
    type (vector2_type), intent(in) :: v1
    real, intent(in) :: fact
    type (vector2_type) :: vec

    vec%x = v1%x * fact
    vec%y = v1%y * fact
  end function 

  function vmag (v1) result(mag)
    type (vector2_type), intent(in) :: v1
    real :: mag

    mag = sqrt((v1%x * v1%x) + (v1%y * v1%y))
  end function 

  function vnormalize (v1) result(vec)
    type (vector2_type), intent(in) :: v1
    real :: mag
    type (vector2_type) :: vec

    mag = sqrt((v1%x**2)+ (v1%y**2))
    vec%x = v1%x / mag
    vec%y = v1%y / mag
  end function 

  function calc_temp (particles, length) result(temp)
    type(particle), dimension(*) :: particles
    integer :: length, i 
    real :: temp

    do i = 1, length
       temp = temp + vmag(particles(i)%speed)
    end do
    temp = sqrt(temp)
  end function calc_temp

  subroutine update_particles (particles, length)
    type(particle), dimension(*) :: particles
    integer :: length
    
    integer :: i
    integer :: ii

    type(vector2_type) :: pos1, pos2
    type (vector2_type) :: speed
    type (vector2_type) :: collision_norm, collision_dir, v1_par, v1_orth, v2_par, v2_orth
    real :: v1_length, v1_length_par, v2_length, v2_length_par, common_vel, v1_length_after, v2_length_after

    do i = 1, length

       pos1 = particles(i)%pos

       speed = particles(i)%speed
 
       do ii = 1, length
          pos2 = particles(ii)%pos

          if (i < ii) then
             if (vmag(vsub(pos1, pos2)) < (2 * MOLEKULA_RADIUS) + 0.1) then

!!$                v1_length = vdot(vsub(particles(i)%speed, particles(ii)%speed), vsub(pos1, pos2))
!!$                v2_length = vdot(vsub(particles(ii)%speed, particles(i)%speed), vsub(pos2, pos1))
!!$                v1_par = vscale(vsub(pos1, pos2), v1_length / (vmag(vsub(pos1, pos2)) ** 2))
!!$                v2_par = vscale(vsub(pos2, pos1),  v2_length / (vmag(vsub(pos2, pos1)) ** 2))
!!$
                ! particles(i)%speed = vscale(vsub(particles(i)%speed, v1_par), ELASTIC_COEFF)
                ! particles(ii)%speed = vscale(vsub(particles(ii)%speed, v2_par), ELASTIC_COEFF)

                collision_norm = vnormalize(vsub(pos1, pos2))
                collision_dir = vector2_type(-collision_norm%y, collision_norm%x)

                v1_length = vdot(particles(i)%speed, collision_norm)
                v2_length = vdot(particles(ii)%speed, collision_norm)

                v1_length_par = vdot(collision_dir, particles(i)%speed)
                v2_length_par = vdot(collision_dir, particles(ii)%speed)

                v1_orth = vscale(collision_norm, (v1_length + v2_length - ELASTIC_COEFF*(v1_length - v2_length)) / 2)
                v2_orth = vscale(collision_norm, (v2_length + v1_length - ELASTIC_COEFF*(v2_length - v1_length)) / 2)

                v1_par = vscale(collision_dir, v1_length_par)
                v2_par = vscale(collision_dir, v2_length_par)

                particles(i)%speed = vadd(v1_par, v1_orth)
                particles(ii)%speed = vadd(v2_par, v2_orth)

                particles(i)%pos = vadd(particles(i)%pos, vscale(particles(i)%speed, dt))
                particles(ii)%pos = vadd(particles(ii)%pos, vscale(particles(ii)%speed, dt))
             end if
          end if

        end do

       if (particles(i)%pos%x >= SCREEN_WIDTH - MOLEKULA_RADIUS .or. particles(i)%pos%x <= MOLEKULA_RADIUS) then
          particles(i)%speed%x = -particles(i)%speed%x
          ! particles(i)%speed = vscale(particles(i)%speed, ELASTIC_COEFF)
       end if

       if (particles(i)%pos%y >= SCREEN_HEIGHT - MOLEKULA_RADIUS .or. particles(i)%pos%y <= MOLEKULA_RADIUS) then
          particles(i)%speed%y = -particles(i)%speed%y
          ! particles(i)%speed = vscale(particles(i)%speed, ELASTIC_COEFF)
       end if

       particles(i)%pos = vadd(particles(i)%pos, vscale(particles(i)%speed, dt))
    end do
 
  end subroutine

  subroutine draw_particle (p)
    type (particle), intent(in) :: p

    integer :: x_pos
    integer :: y_pos

    x_pos = int(p%pos%x)
    y_pos = int(p%pos%y)

    if (x_pos < SCREEN_WIDTH .and. y_pos < SCREEN_HEIGHT) then
       call draw_circle(x_pos, y_pos, MOLEKULA_RADIUS, BLUE)
    end if
  end subroutine

  subroutine draw_particles (particles, n)
    type(particle), dimension(*), intent(in) :: particles
    integer, intent(in) :: n

    integer :: i

    do i = 1, n
       call draw_particle(particles(i))
    end do

  end subroutine

  subroutine print_particles (particles, length)
    type(particle), dimension(*) :: particles
    integer :: length
    integer :: i

   end subroutine

  
end program
   
