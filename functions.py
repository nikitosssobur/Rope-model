import json
import matplotlib.pyplot as plt
import numpy as np
from prettytable import PrettyTable


def create_model_config(config_data):
    output_config = dict()
    output_config['x_coords'], output_config['y_coords'] = config_data['x_coords'], config_data['y_coords']
    output_config['radius'] = config_data['fixed_point_radius']
    output_config['moving_points_num'] = config_data['hor_rope_moving_points_num']
    output_config['moving_points_radius'] = config_data['hor_rope_moving_points_radius']
    output_config['moving_points_velocities'] = config_data['hor_rope_moving_points_velocities']
    output_config['rope_mass'] = config_data['hor_rope_mass']
    output_config['tension_force'] = config_data['tension_force']
    output_config['stiffness'] = config_data['hor_rope_stiffness']
    output_config['viscosity'] = config_data['hor_rope_viscosity']
    output_config['rest_len'] = config_data['ver_rope_rest_len']
    output_config['point_num'] = config_data['ver_rope_points_num']
    output_config['weighted_point_mass'] = config_data['weighted_point_mass']
    output_config['weighted_point_radius'] = config_data['weighted_point_radius']
    output_config['ver_rope_point_radius'] = config_data['ver_rope_point_radius']
    output_config['ver_rope_mass'] = config_data['ver_rope_mass']
    output_config['ver_rope_stiffness'] = config_data['ver_rope_stiffness']
    output_config['ver_rope_viscosity'] = config_data['ver_rope_viscosity']
    output_config['attach_point_id'] = config_data['attach_point_id']
    output_config['surface_color'] = config_data['surface_color']
    output_config['display_width'] = config_data['display_width']
    output_config['display_height'] = config_data['display_height']
    output_config['print_system_info'] = config_data['print_system_info']
    output_config['real_meter_dist'] = config_data['real_meter_dist']
    output_config['fall_height'] = config_data['fall_height']
    output_config['sec_num'] = config_data['sec_num']
    output_config['time_steps_per_sec'] = config_data['time_steps_per_sec']
    output_config['solver_method'] = config_data['solver_method']
    output_config['rope_fragments_indexes'] = config_data['rope_fragments_indexes']
    return output_config


def get_info(valid_args, rope):
    fixed_points_table = PrettyTable()
    fixed_points_table.title = 'Fixed points data'
    fixed_points_table.field_names = ["x coord (in pixels)", "y coord (in pixels)", "x coord (in meters)", "y coord (in meters)", "radius (in pixels)"]
    fixed_points_table.add_row([valid_args['x_coords'][0], valid_args['y_coords'][0], 
                                rope.objects['fixed_points'][0].si_point_data['coordinates'][0], 
                                rope.objects['fixed_points'][0].si_point_data['coordinates'][1], 
                                valid_args['radius']])
        
    fixed_points_table.add_row([valid_args['x_coords'][1], valid_args['y_coords'][1], 
                                rope.objects['fixed_points'][1].si_point_data['coordinates'][0], 
                                rope.objects['fixed_points'][1].si_point_data['coordinates'][1], 
                                valid_args['radius']])
        
    ropes_table = PrettyTable()
    ropes_table.title = 'Ropes data'
    ropes_table.add_column(fieldname="Info", column=["points number", "rope mass (in kg)", "tension force (Newtons)", 
                                "stiffness (Newtons)", "viscosity (Newtons forces)", "rest lenght (in pixels)", "rest lenght (in meters)", 
                                "moving point radius (in pixels)", "moving point mass (in kg)", "attachment point index"])
    
    ropes_table.add_column(fieldname="Horizontal", column=[valid_args['moving_points_num'], valid_args['rope_mass'], 
                               valid_args['tension_force'], valid_args['stiffness'], valid_args['viscosity'], 
                               rope.hor_rest_len*rope.pix_per_metr, rope.hor_rest_len, 
                               valid_args['moving_points_radius'], valid_args['rope_mass']/valid_args['moving_points_num'], 
                               valid_args['attach_point_id']]) 
    
    ropes_table.add_column(fieldname="Vertical", column=[valid_args['point_num'], valid_args['ver_rope_mass'], 0, 
                               valid_args['ver_rope_stiffness'], valid_args['ver_rope_viscosity'], 
                               valid_args['rest_len']*rope.pix_per_metr, valid_args['rest_len'], 
                               valid_args['ver_rope_point_radius'], valid_args['ver_rope_mass']/valid_args['point_num'], 
                               valid_args['attach_point_id']])
        

    window_table = PrettyTable()
    window_table.title = 'Window data'
    window_table.field_names = ["Display width", "Display height", 'Surface color']
    window_table.add_row([valid_args['display_width'], valid_args['display_height'], valid_args['surface_color']])
        
 
    experiment_table = PrettyTable()
    experiment_table.title = 'Experiment data'
    experiment_table.add_column(fieldname="Info", column=["Distance between fixed points (in meters)",
                                "Weighted point fall height (in meters)", "Pixels per meter coefficient", 
                                "Experiment duration (in seconds)", "Time steps per second",
                                "ODE Solver method", "Rope fragments indexes (for force graph)"])
    
    experiment_table.add_column(fieldname="Values", column=[valid_args['real_meter_dist'], 
                                valid_args['fall_height'], rope.pix_per_metr, valid_args['sec_num'], 
                                valid_args['time_steps_per_sec'], valid_args['solver_method'], 
                                [valid_args['rope_fragments_indexes'] if valid_args['rope_fragments_indexes'] is not None 
                                 else valid_args['moving_points_num']]])

    return fixed_points_table, ropes_table, window_table, experiment_table 


def unpack_model_config(config):
    if isinstance(config, str) and config[-5:] == '.json':
        with open(config, 'r') as f:
            loaded_config = json.load(f)
    elif isinstance(config, dict):
        loaded_config = config 
    else:
        raise NotImplementedError

    fixed_point_data = {'x_coords': loaded_config['x_coords'], 'y_coords': loaded_config['y_coords'], 
                        'radius': loaded_config['radius']}
    
    moving_point_data = {'moving_points_num': loaded_config['moving_points_num'], 
                         'moving_points_radius': loaded_config['moving_points_radius'], 
                         'moving_points_velocities': loaded_config['moving_points_velocities']}
    
    hor_rope_data = {'fixed_point_data': fixed_point_data, 'moving_point_data': moving_point_data,
                    'rope_mass': loaded_config['rope_mass'], 'rope_stiffness': loaded_config['stiffness'], 
                    'rope_viscosity': loaded_config['viscosity'], 
                    'tension_force': loaded_config['tension_force']}
    
    ver_rope_data = {'rest_len': loaded_config['rest_len'], 
                'point_num': loaded_config['point_num'], 
                'weighted_point_mass': loaded_config['weighted_point_mass'], 
                'weighted_point_radius': loaded_config['weighted_point_radius'],
                'ver_rope_point_radius': loaded_config['ver_rope_point_radius'],
                'rope_mass': loaded_config['ver_rope_mass'],
                'rope_stiffness': loaded_config['ver_rope_stiffness'],
                'rope_viscosity': loaded_config['ver_rope_viscosity']}
    
    surface_data = {'color': tuple(loaded_config['surface_color']), 
                    "width": loaded_config['display_width'], 
                    "height": loaded_config['display_height']}
    
    additional_info = {'attach_point_id': loaded_config['attach_point_id'], 
                       'print_system_info': loaded_config['print_system_info'], 
                       'real_meter_dist': loaded_config['real_meter_dist'], 
                       'fall_height': loaded_config['fall_height'], 
                       'sec_num': loaded_config['sec_num'], 
                       'time_steps_per_sec': loaded_config['time_steps_per_sec'], 
                       'method': loaded_config['solver_method'], 
                       'rope_fragments_indexes': loaded_config['rope_fragments_indexes']}
    
    return hor_rope_data, ver_rope_data, surface_data, additional_info


def draw_force_graph(force_vals, time_steps=None, y_label='F(t)', x_label='t'):
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    if time_steps is None:
        plt.plot(force_vals)
    else:
        plt.plot(time_steps, force_vals)
    plt.show()


def find_time_steps(sec_num=60, time_steps_per_second=100):
    return {'general_time_steps': sec_num * time_steps_per_second, 'delay': int(1000 / time_steps_per_second)}




