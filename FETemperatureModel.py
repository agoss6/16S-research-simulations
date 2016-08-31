__author__ = 'Sandy'

# This python file simulates the temperature of square copper wire with current and water running through it
# as a function of time. The parameters that can be tweaked are: the length of the wire, the current,
# the temperature of the water that enters the wire, the velocity of the water, and the time resolution of
# the simulation.
#
# The idea is that the copper tube can be thought of as many smaller elements of length dx.
# The relationship dx = flow_velocity * dt  must hold so dx is calculated from that. Each element has 1.energy
# added to the copper due to Ohmic heating, 2.energy lost to the water flowing through it, and 3.that same
# energy added to the water flowing through it. 1 can be calculating using power = I^2*R. 2 and 3 can be calculated by
# assuming that each element has a uniform heat flux, allowing us to use the formula
# (heat transfer rate) = (coefficient of heat transfer)*(area contacting the water)*(temperature difference)
#
# The code works by initializing two python lists that represent the temperature gradient of the copper wire
# and the water in the copper wire. Each time step dt, a while loop goes through each element in the gradients.
# It adds the energy due to ohmic heating to the copper, then transfers the energy from the copper to the water.
# Finally, it removes the last element in the water gradient and adds a new element to the front of the water gradient
# to simulate water flowing through the tube.

# properties of the copper tube
copper_density = 8915.0 # kilograms/meter^3  (average density of the copper = (8890+8940)/2 )
copper_cp = 385.0 # Joules / (kilogram * Kelvin)  (specific heat of copper at constant pressure)

########VARIABLE#############
length = 1.0 # meter (length of wire)
#####################
outer_side = .00375  # meter (1/8 inch, the outer side length)
inner_bore = .00127  # meter (1/20 inch, the inner bore side length)
tube_area = outer_side**2 - inner_bore**2 # meter^2  (cross sectional area of the square tubing)

alphaTemp = 3.86e-3 # K^-1  (temperature coefficient of resistance for copper)
rho0 = 1.71e-8  # ohm * meter  (resistivity at room temperature)
r0 = rho0 / tube_area  # ohms/meter  (resistance per meter at room temperature)



# properties of the water
########VARIABLE#############
water_inlet_temperature = 323.0
#####################
water_density = 1000.0 # kg/m^3
water_cp = 4179.0 # Joules/ (kilogram*Kelvin) (specific heat of water at constant pressure)
########VARIABLE#############
flow_velocity = 1.0  # m/s
#####################


# sizes of each finite element
########VARIABLE#############
dt = .01 # time resolution (seconds)
#####################
dx = flow_velocity * dt  # length of each tube segment (dependent on velocity and time reolution)


# other constants
room_temperature = 293 # Kelvin  (room temperature)
########VARIABLE#############
I = 500.0 # Amperes (current through the wire)
#####################


# This is a function that calculates the resistance per meter of square copper tubing as a function of temperature.
# It takes into account how resistivity changes with temperature using the temperature coefficient.
# equation used: R(T) = R0 [1 + alphaTemp(T- T0)]
def resistance(T):
    unit_resistance = (r0*(1+alphaTemp*(T-room_temperature)))
    return unit_resistance



# Simulation of the temperature as a function of time.
def temperature(time):
    if time < 0:
        return "Time cannot be negative."

    #initialize copper and water gradients
    starting_temperature = water_inlet_temperature
    copperGradient = []
    waterGradient = []
    gradient_size = int(length / dx)
    for i in range(gradient_size):
        copperGradient.append(starting_temperature)
        waterGradient.append(starting_temperature)

    time_elapsed = 0.0 # time elapsed

    # while loop that runs for the specified time
    while time_elapsed < time:
        i = 0
        # while loop that goes through each element in the gradients and calculates the energy gains and losses to each one
        while i < gradient_size:
            # energy gained from ohmic heating
            current_dE = (I*I) * resistance(copperGradient[i])* dx * dt  # change in internal energy due to ohmic heating
            current_dT = current_dE / (copper_cp * dx * tube_area * copper_density)  # converting energy change to temperature change
            copperGradient[i] += current_dT

            # heat transfer rate between the copper and the water in an element
            h = 3000 # watts/kelvin (convective heat transfer coefficient for copper and water given by engineering toolbox)
            A = 4 * dx * inner_bore # m^2 (contact area between the copper and the water)
            deltaT = copperGradient[i] - waterGradient[i] # temperature difference
            heat_transfer_rate = h*A*deltaT # joules/second (convective heat transfer rate using q = h * A * deltaT)

            # energy and temperature changes due to convective cooling
            convective_dE = heat_transfer_rate * dt
            copper_dT = -convective_dE / (copper_cp * dx * tube_area * copper_density)
            water_dT = convective_dE / (water_cp * dx * tube_area * water_density)
            copperGradient[i] += copper_dT
            waterGradient[i] += water_dT

            i += 1

        # simulate the water flowing
        waterGradient.insert(0, water_inlet_temperature)
        waterGradient.pop()
        time_elapsed += dt

    return copperGradient
    # return waterGradient
