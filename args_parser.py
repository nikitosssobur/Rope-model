import argparse as ap
import json 
from model_config_creator import DEFAULT_CONFIG_PATH


with open(DEFAULT_CONFIG_PATH, "r") as f:
    default_config = json.load(f)



parser = ap.ArgumentParser(description="Simulation model")

parser.add_argument('-xs', '--x_coords', type=float, nargs='+', 
                    default=default_config['x_coords'], 
                    help="x coords of the fixed left and right edge points in the window (in pixels)")

parser.add_argument('-ys', '--y_coords', type=float, nargs='+', 
                    default=default_config['y_coords'], 
                    help="y coords of the fixed left and right edge points in the window (in pixels)")

parser.add_argument('-fpr', '--fixed_point_radius', type=float, 
                    default=default_config['radius'], 
                    help="radius of the fixed edge points (in pixels)")

parser.add_argument('-hrmpn', '--hor_rope_moving_points_num', type=int, 
                    default=default_config['moving_points_num'],
                    help="horizontal rope moving points number (integer, including ropes attachment point)")

parser.add_argument('-hrmpr', '--hor_rope_moving_points_radius', type=float,
                    default=default_config['moving_points_radius'],
                    help="horizontal rope moving points radius (in pixels)")

parser.add_argument('-hrm', '--hor_rope_mass', type=float,
                    default=default_config['rope_mass'], 
                    help="total mass of the whole horizontal rope (in kg, mass of all points of the rope)")

parser.add_argument('-hrmpv', '--hor_rope_moving_points_velocities', type=float, nargs='+', 
                    default=default_config['moving_points_velocities'], 
                    help="initial velocity vector of all points of the horizontal rope (default: [0, 0], m/sec)")

parser.add_argument('-tf', '--tension_force', type=float, default=default_config['tension_force'],
                    help="initial tension force in horizontal rope (Newtons)")

parser.add_argument('-hrs', '--hor_rope_stiffness', type=float, default=default_config['stiffness'],
                    help="stiffness coefficient of the horizontal rope")

parser.add_argument('-hrv', '--hor_rope_viscosity', type=float, default=default_config['viscosity'],
                    help="viscosity coefficient of the horizontal rope")

parser.add_argument('-vrrl', '--ver_rope_rest_len', type=float, default=default_config['rest_len'],
                    help="vertical rope rest length (in meters)")

parser.add_argument('-vrpn', '--ver_rope_points_num', type=int, default=default_config['point_num'],
                    help="vertical rope points number (attachment point excluded, weighted point included)")

parser.add_argument('-wpm', '--weighted_point_mass', type=float, default=default_config['weighted_point_mass'],
                    help="weighted point mass (in kg)")

parser.add_argument('-wpr', '--weighted_point_radius', type=float, default=default_config['weighted_point_radius'],
                    help="weighted point radius (in pixels)")

parser.add_argument('-vrpr', '--ver_rope_point_radius', type=float, default=default_config['ver_rope_point_radius'],
                    help="vertical rope point radius (in pixels)")

parser.add_argument('-vrm', '--ver_rope_mass', type=float,
                    default=default_config['ver_rope_mass'], 
                    help="total mass of the whole vertical rope (in kg, mass of all points of the rope)")

parser.add_argument('-vrs', '--ver_rope_stiffness', type=float, default=default_config['ver_rope_stiffness'],
                    help="stiffness coefficient of the vertical rope")

parser.add_argument('-vrv', '--ver_rope_viscosity', type=float, default=default_config['ver_rope_viscosity'],
                    help="viscosity coefficient of the vertical rope")

parser.add_argument('-apid', '--attach_point_id', type=int, default=default_config['attach_point_id'], 
                    help="index of the attachment moving point of horizontal and vertical ropes (points are numbered from 0, from the left to the right)")

parser.add_argument('-sc', '--surface_color', type=int, nargs='+', default=default_config['surface_color'],
                    help="surface color in RGB format (default white: [255, 255, 255])")

parser.add_argument('-dw', '--display_width', type=int, default=default_config['display_width'], 
                    help="display width (integer in pixels)")

parser.add_argument('-dh', '--display_height', type=int, default=default_config['display_height'],
                    help="display height (integer in pixels)")

parser.add_argument('-psi', '--print_system_info', action='store_false', default=default_config['print_system_info'], 
                    help="flag (bool) for printing info about the rope system in console (if argument passed, info will not be printed (False), else info printed, default: True)")

parser.add_argument('-rmd', '--real_meter_dist', type=float, default=default_config['real_meter_dist'], 
                    help="real meter distance between fixed points (default: 100 meters)")

parser.add_argument('-fh', '--fall_height', type=float, default=default_config['fall_height'], 
                    help="the height of the point from which the weighted point falls (meters, used for initial velocities of vertical rope points for sqrt(2*g*h))")

parser.add_argument('-sn', '--sec_num', type=float, default=default_config['sec_num'], 
                    help="duration of the rope simulation (in seconds)")

parser.add_argument('-tsps', '--time_steps_per_sec', type=int, default=default_config['time_steps_per_sec'], 
                    help="number of time steps per second (fps)")

parser.add_argument('-sm', '--solver_method', type=str, default=default_config['solver_method'], 
                    help='''numerical solution method of the ODE system (namings of available methods from scipy, default: Radau), " \
                    also available RK45, RK23, etc. Additional info about methods can be found here: 
                    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#r179348322575-4''')

parser.add_argument('-rfi', '--rope_fragments_indexes', type=int, nargs='+', 
                    default=default_config['rope_fragments_indexes'], 
                    help="list of indexes of the ropes for which forces graphs will be builded (by default only for right edge rope fragment connected between last horizonal rope moving point and right fixed point graph is builded)")