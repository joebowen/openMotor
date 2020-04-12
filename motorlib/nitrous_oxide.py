# Tank emptying code extracts
# (c) Rick Newlands, AspireSpace 
# Variables are in metric units except where stated otherwise 

import math

pcrit = 72.51  # critical pressure, Bar Abs
rho_crit = 452.0  # critical density, kg/m3
tcrit = 309.57  # critical temperature, Kelvin (36.42 Centigrade)
zcrit = 0.28  # critical compressibility factor
gamma = 1.3  # average over subcritical range

CENTIGRADE_TO_KELVIN = 273.15  # to Kelvin

# ranges of function validity
lower_temp_limit = -90.0 + CENTIGRADE_TO_KELVIN
upper_temp_limit = 36.4 + CENTIGRADE_TO_KELVIN

BAR_TO_PASCALS = 100000.0
PASCALS_TO_BAR = (1.0 / BAR_TO_PASCALS)


def sqr(bob):
    return bob * bob


def sgn(bob):
    # signum of a number, used below
    if bob >= 0.0:
        signum = 1
    else:
        signum = -1

    return signum


# Nitrous oxide vapour pressure, Bar
def nox_vp(t_kelvin):
    p = (
        1.0,
        1.5,
        2.5,
        5.0
    )
    b = (
        -6.71893,
        1.35966,
        -1.3779,
        -4.051
    )

    tr = t_kelvin / tcrit
    rab = 1.0 - tr
    shona = 0.0

    for dd in range(0, 4):
        shona += b[dd] * math.pow(rab, p[dd])

    bob = pcrit * math.exp((shona / tr))

    return bob


# Nitrous oxide saturated liquid density, kg/m3
def nox_lrho(t_kelvin):
    b = (
        1.72328,
        -0.8395,
        0.5106,
        -0.10412
    )

    tr = t_kelvin / tcrit
    rab = 1.0 - tr
    shona = 0.0

    for dd in range(0, 4):
        shona += b[dd] * math.pow(rab, ((dd + 1) / 3.0))

    bob = rho_crit * math.exp(shona)

    return bob


# Nitrous oxide saturated vapour density, kg/m3
def nox_vrho(t_kelvin):
    b = (
        -1.009,
        -6.28792,
        7.50332,
        -7.90463,
        0.629427
    )

    tr = t_kelvin / tcrit
    rab = (1.0 / tr) - 1.0
    shona = 0.0

    for dd in range(0, 5):
        shona += b[dd] * math.pow(rab, ((dd + 1) / 3.0))

    bob = rho_crit * math.exp(shona)

    return bob


# Nitrous liquid Enthalpy (Latent heat) of vaporisation, J/kg
def nox_enth_v(t_kelvin):
    b_l = (
        -200.0,
        116.043,
        -917.225,
        794.779,
        -589.587
    )

    b_v = (
        -200.0,
        440.055,
        -459.701,
        434.081,
        -485.338
    )

    tr = t_kelvin / tcrit
    rab = 1.0 - tr
    shona_l = b_l[0]
    shona_v = b_v[0]

    for dd in range(1, 5):
        shona_l += b_l[dd] * math.pow(rab, (dd / 3.0))  # saturated liquid enthalpy
        shona_v += b_v[dd] * math.pow(rab, (dd / 3.0))  # saturated vapour enthalpy

    bob = (shona_v - shona_l) * 1000.0  # net during change from liquid to vapour

    return bob


# Nitrous saturated liquid isobaric heat capacity, J/kg K
def nox_cpl(t_kelvin):
    b = (
        2.49973,
        0.023454,
        -3.80136,
        13.0945,
        -14.518
    )

    tr = t_kelvin / tcrit
    rab = 1.0 - tr
    shona = 1.0 + b[1] / rab

    for dd in range(1, 4):
        shona += b[(dd + 1)] * math.pow(rab, dd)

    bob = b[0] * shona * 1000.0  # convert from KJ to J

    return bob


# liquid nitrous thermal conductivity, W/m K
def nox_kl(t_kelvin):
    b = (
        72.35,
        1.5,
        -3.5,
        4.5
    )

    # max. 10 deg C
    if t_kelvin > 283.15:
        tr = 283.15 / tcrit
    else:
        tr = t_kelvin / tcrit

    rab = 1.0 - tr
    shona = 1.0 + b[3] * rab

    for dd in range(1, 3):
        shona += b[dd] * math.pow(rab, (dd / 3.0))

    bob = b[0] * shona / 1000  # convert from mW to W

    return bob


# nitrous temperature based on pressure (bar)
def nox_on_press(p_bar_abs):
    p = (
        1.0,
        1.5,
        2.5,
        5.0
    )

    b = (
        -6.71893,
        1.35966,
        -1.3779,
        -4.051
    )

    step = -1.0
    temp_k = (tcrit - 0.1) - step

    # iterative loop
    while True:
        while True:
            temp_k += step
            tr = temp_k / tcrit
            rab = 1.0 - tr
            shona = 0.0

            for dd in range(0, 4):
                shona += b[dd] * math.pow(rab, p[dd])

            pp_guess = pcrit * math.exp(shona / tr)

            if not ((pp_guess - p_bar_abs) * sgn(step)) < 0.0:
                break

        if not math.fabs((pp_guess - p_bar_abs)) > 0.01:
            break

        step = step / (-2.0)  # reduce step size

    bob = temp_k

    return bob  # return temperature


def linear_interpolate(x, x1, y1, x2, y2):
    # Linear interpolation routine, with limiters added
    # incase x isn't within the range range x1 to x2
    # Limits updated to allow descending x values
    # x = input value
    # x1 = MINIMUM BOUNDS ( minimum x )
    # y1 = MINIMUM VALUE ( output value at MINIMUM BOUNDS )
    # x2 = MAXIMUM BOUNDS ( maximum x )
    # y2 = MAXIMUM VALUE ( output value at MAXIMUM BOUNDS )
    ##

    # This procedure extrapolates the y value for the x position
    # on a line defined by x1,y1; x2,y2
    # the constants to find in y=mx+b

    if (x1 < x2) and ((x <= x1) or (x >= x2)):  # ascending x values
        if x <= x1:
            return y1
        else:
            return y2
    elif (x1 > x2) and ((x >= x1) or (x <= x2)):  # descending x values
        if x >= x1:
            return y1
        else:
            return y2
    else:
        m = (y2 - y1) / (x2 - x1)  # calculate the gradient m
        c = y1 - m * x1  # calculate the y-intercept c
        y = m * x + c  # the final calculation

    return y


def compressibility_factor(p_bar_abs):
    # Compressibility factor of vapour (subcritical) on the saturation line
    return linear_interpolate(p_bar_abs, 0.0, 1.0, pcrit, zcrit)


class N2OTank:
    def __init__(self, delta_time, initial_tank_pressure, initial_ullage_percentage, tank_volume, num_of_injectors, injector_diameter, injector_k2_coefficient):
        # Subroutine to initialise main program variables
        # Gets called only once, prior to firing

        self.nox_pcrit = 72.51
        self.nox_zcrit = 0.28
        self.nox_gamma = 1.3

        self.tank_volume = tank_volume

        self.delta_time = delta_time

        self.vapor_phase = False
        self.hybrid_fault = 0
        self.tank_vapour_mass = 0.0
        self.mdot_tank_outflow = 0.0
        self.old_mdot_tank_outflow = 0.0

        # set initial nitrous vapour (tank) pressure
        self.tank_fluid_temperature_K = nox_on_press(initial_tank_pressure)  # set tank pressure
        self.tank_pressure = initial_tank_pressure

        # reality check
        if self.tank_fluid_temperature_K > (36.0 + CENTIGRADE_TO_KELVIN):
            self.tank_fluid_temperature_K = 36.0 + CENTIGRADE_TO_KELVIN
            self.hybrid_fault = 2

        # get initial nitrous properties
        self.tank_liquid_density = nox_lrho(self.tank_fluid_temperature_K)
        self.tank_vapour_density = nox_vrho(self.tank_fluid_temperature_K)

        # base the nitrous vapour volume on the tank percentage ullage (gas head-space)
        self.tank_vapour_volume = (initial_ullage_percentage / 100.0) * tank_volume
        self.tank_liquid_volume = tank_volume - self.tank_vapour_volume
        self.tank_liquid_mass = self.tank_liquid_density * self.tank_liquid_volume
        self.tank_vapour_mass = self.tank_vapour_density * self.tank_vapour_volume
        self.tank_propellant_contents_mass = self.tank_liquid_mass + self.tank_vapour_mass  # total mass within tank
        self.prev_tank_propellant_contents_mass = self.tank_propellant_contents_mass

        # initialise values needed later
        self.old_liquid_nox_mass = self.tank_liquid_mass
        self.old_vapour_nox_mass = self.tank_vapour_mass
        self.initial_liquid_propellant_mass = self.tank_liquid_mass
        self.initial_vapour_propellant_mass = self.tank_vapour_mass

        # guessed initial value of amount of nitrous vaporised per iteration
        # in the nitrous tank blowdown model (actual value is not important)
        self.vaporised_mass_old = 0.001

        # individual injector orifice total loss coefficent K2
        bob = math.pi * sqr((injector_diameter / 2.0))  # orifice cross sectional area

        self.injector_loss_coefficient = (injector_k2_coefficient / (sqr((num_of_injectors * bob)))) * PASCALS_TO_BAR

        self.initial_vapour_temp_k = self.tank_fluid_temperature_K
        self.initial_vapour_mass = self.tank_vapour_mass
        self.initial_vapour_pressure_bar = initial_tank_pressure
        self.initial_vapour_density = self.tank_vapour_density
        self.initial_z = compressibility_factor(initial_tank_pressure)

        self.vapour_temp_k = self.tank_fluid_temperature_K
        self.z = compressibility_factor(initial_tank_pressure)

    def injector_model(self, upstream_pressure, downstream_pressure):
        # calculate injector pressure drop (Bar) and mass flowrate (kg/sec)

        pressure_drop = upstream_pressure - downstream_pressure  # Bar

        # reality check
        if pressure_drop < 0.00001:
            pressure_drop = 0.00001

        # is injector pressure drop lower than 20 percent of chamber pressure?
        if (pressure_drop / downstream_pressure) < 0.2:
            self.hybrid_fault = 3  # too low for safety

        # Calculate fluid flowrate through the injector, based on the
        # total-pressure loss factor between the tank and combustion chamber
        # (injector_loss_coefficient includes K coefficient and orifice cross-sectional areas)
        mass_flowrate = math.sqrt((2.0 * self.tank_liquid_density * pressure_drop / self.injector_loss_coefficient))

        return mass_flowrate  # kg/sec

    def get_mass_flow_per_iteration(self, chamber_pressure_bar):
        if self.hybrid_fault == 9:
            self.vapor_phase = True

            self.initial_vapour_temp_k = self.tank_fluid_temperature_K
            self.initial_vapour_mass = self.tank_vapour_mass
            self.initial_vapour_pressure_bar = self.tank_pressure
            self.initial_vapour_density = self.tank_vapour_density
            self.initial_z = compressibility_factor(self.tank_pressure)

            self.vapour_temp_k = self.tank_fluid_temperature_K
            self.z = compressibility_factor(self.tank_pressure)

        if not self.vapor_phase:
            self.nitrous_tank_liquid(chamber_pressure_bar)
        else:
            self.subcritical_tank_no_liquid(chamber_pressure_bar)

        mass_flows = self.prev_tank_propellant_contents_mass - self.tank_propellant_contents_mass

        self.prev_tank_propellant_contents_mass = self.tank_propellant_contents_mass

        if self.tank_propellant_contents_mass < 0:
            print('self.tank_propellant_contents_mass < 0')
            mass_flows = 0

        print(f'tank_pressure (bar): {self.tank_pressure}')

        if self.tank_pressure < chamber_pressure_bar:
            print('self.tank_pressure < chamber_pressure_bar')
            mass_flows = 0

        return mass_flows

    def nitrous_tank_liquid(self, chamber_pressure_bar):
        # Equilibrium (instantaneous boiling) tank blowdown model
        # Empty tank of liquid nitrous

        lagged_bob = 0.0

        # blowdown simulation using nitrous oxide property calcs subroutines
        # update last-times values, O = 'old'
        omdot_tank_outflow = self.mdot_tank_outflow
        enth_of_vap = nox_enth_v(self.tank_fluid_temperature_K)

        # Get enthalpy (latent heat) of vaporisation
        spec_heat_cap = nox_cpl(self.tank_fluid_temperature_K)

        # Get specific heat capacity of the liquid nitrous
        # Calculate the heat removed from the liquid nitrous during its vaporisation
        deltaq = self.vaporised_mass_old * enth_of_vap

        # temperature drop of the remaining liquid nitrous due to losing this heat
        deltatemp = -(deltaq / (self.tank_liquid_mass * spec_heat_cap))
        self.tank_fluid_temperature_K += deltatemp  # update fluid temperature

        # reality checks
        if self.tank_fluid_temperature_K < (-90.0 + CENTIGRADE_TO_KELVIN):
            self.tank_fluid_temperature_K = (-90.0 + CENTIGRADE_TO_KELVIN)  # lower limit
            self.hybrid_fault = 1
        elif self.tank_fluid_temperature_K > (36.0 + CENTIGRADE_TO_KELVIN):
            self.tank_fluid_temperature_K = (36.0 + CENTIGRADE_TO_KELVIN)  # upper limit
            self.hybrid_fault = 2

        # get current nitrous properties
        self.tank_liquid_density = nox_lrho(self.tank_fluid_temperature_K)
        self.tank_vapour_density = nox_vrho(self.tank_fluid_temperature_K)
        self.tank_pressure = nox_vp(self.tank_fluid_temperature_K)

        # calculate injector pressure drop and mass flowrate
        mdot_tank_outflow = self.injector_model(self.tank_pressure, chamber_pressure_bar)

        # integrate mass flowrate using Addams second order integration formula
        # (my preferred integration formulae, feel free to choose your own.)
        # Xn=X(n-1) + DT/2 * ((3 * Xdot(n-1) - Xdot(n-2))
        # O infront of a variable name means value from previous timestep (Old)
        delta_outflow_mass = 0.5 * self.delta_time * (3.0 * mdot_tank_outflow - omdot_tank_outflow)

        # drain the tank based on flowrates only
        self.tank_propellant_contents_mass -= delta_outflow_mass

        # update mass within tank for next iteration
        self.old_liquid_nox_mass -= delta_outflow_mass

        # update liquid mass within tank for next iteration

        # now the additional effects of phase changes
        # The following equation is applicable to the nitrous tank, containing saturated nitrous:
        # tank_volume = liquid_nox_mass / liquid_nox_density + gaseous_nox_mass / gaseous_nox_density
        # Rearrange this equation to calculate current liquid_nox_mass
        bob = (1.0 / self.tank_liquid_density) - (1.0 / self.tank_vapour_density)

        self.tank_liquid_mass = (self.tank_volume - (self.tank_propellant_contents_mass / self.tank_vapour_density)) / bob
        self.tank_vapour_mass = self.tank_propellant_contents_mass - self.tank_liquid_mass

        # update for next iteration
        bob = self.old_liquid_nox_mass - self.tank_liquid_mass

        # Add a 1st-order time lag (of 0.15 seconds) to aid numerical
        # stability (this models the finite time required for boiling)
        tc = self.delta_time / 0.15
        lagged_bob = tc * (bob - lagged_bob) + lagged_bob  # 1st-order lag
        self.vaporised_mass_old = lagged_bob

        # Check for model fault at nearly zero liquid oxidiser mass
        # If this occurs, use the fault flag to trigger burnout
        if self.tank_liquid_mass > self.old_liquid_nox_mass:
            self.hybrid_fault = 9

        # update tank contents for next iteration
        self.old_liquid_nox_mass = self.tank_liquid_mass

    def subcritical_tank_no_liquid(self, chamber_pressure_bar):
        # subroutine to model the tank emptying of vapour only
        # Isentropic vapour-only blowdown model

        # calculate injector pressure drop and mass flowrate
        self.mdot_tank_outflow = self.injector_model(self.tank_pressure, chamber_pressure_bar)

        # integrate mass flowrate using Addams second order integration formula
        # Xn=X(n-1) + DT/2 * ((3 * Xdot(n-1) - Xdot(n-2))
        delta_outflow_mass = 0.5 * self.delta_time * (3.0 * self.mdot_tank_outflow - self.old_mdot_tank_outflow)

        # drain the tank based on flowrates only
        self.tank_propellant_contents_mass -= delta_outflow_mass  # update mass within tank for next iteration

        # drain off vapour
        self.tank_vapour_mass -= delta_outflow_mass  # update vapour mass within tank for next iteration

        # initial guess
        current_z_guess = compressibility_factor(self.tank_pressure)
        step = 1.0 / 0.9  # initial step size
        aim = 0  # flags used below to home-in

        # recursive loop to get correct compressibility factor
        while True:
            # use isentropic relationships
            bob = self.nox_gamma - 1.0
            self.vapour_temp_k = self.initial_vapour_temp_k * math.pow(((self.tank_vapour_mass * current_z_guess) / (self.initial_vapour_mass * self.initial_z)), bob)

            bob = self.nox_gamma / (self.nox_gamma - 1.0)
            self.tank_pressure = self.initial_vapour_pressure_bar * math.pow((self.vapour_temp_k / self.initial_vapour_temp_k), bob)

            self.z = compressibility_factor(self.tank_pressure)

            old_aim = aim

            if current_z_guess < self.z:
                current_z_guess *= step
                aim = 1

            else:
                current_z_guess /= step
                aim = -1

            # check for overshoot of target, and if so, reduce step nearer to 1.0
            if aim == -old_aim:
                step = math.sqrt(step)

            if not (((current_z_guess / self.z) > 1.000001) or ((current_z_guess / self.z) < (1.0 / 1.000001))):
                break

        bob = 1.0 / (self.nox_gamma - 1.0)

        self.tank_vapour_density = self.initial_vapour_density * math.pow((self.vapour_temp_k / self.initial_vapour_temp_k), bob)


if __name__ == '__main__':
    delta_time = 0.1
    initial_tank_pressure = 49
    initial_ullage_percentage = 5
    tank_volume = 100
    num_of_injectors = 1
    injector_diameter = 0.1
    injector_k2_coefficient = 0.1

    tank = N2OTank(
        delta_time,
        initial_tank_pressure,
        initial_ullage_percentage,
        tank_volume,
        num_of_injectors,
        injector_diameter,
        injector_k2_coefficient
    )

    chamber_pressure = 34

    mass_flow = tank.get_mass_flow_per_iteration(chamber_pressure)

    while mass_flow > 0:
        mass_flow = tank.get_mass_flow_per_iteration(chamber_pressure)

        print(f'{mass_flow}, {tank.tank_propellant_contents_mass}, {tank.tank_pressure} {tank.hybrid_fault} {tank.vapor_phase}')
