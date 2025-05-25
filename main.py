import pygame as pg
from data_validation import ArgsValidator
from args_parser import parser
from functions import create_model_config, unpack_model_config, find_time_steps, get_info
from simulation_class import Simulation
from rope_model import RopeModel



if __name__=='__main__':
    parsed_args = parser.parse_args()
    temp_config = create_model_config(vars(parsed_args))
    valid_args = dict(ArgsValidator(**temp_config))
    print('Input arguments successfully validated!')
    
    
    hor_rope_data, ver_rope_data, surface_data, additional_info = unpack_model_config(valid_args) 

    pg.init()
    surface = pg.display.set_mode((surface_data['width'],  surface_data['height']))

    rope = RopeModel(surface, hor_rope_data=hor_rope_data, real_meter_dist=additional_info['real_meter_dist'], 
                    fall_height=additional_info['fall_height'])
    rope.set_ver_rope_data(ver_rope_data)
    rope.create_model()
    rope.create_free_end_rope(attach_point_id=additional_info['attach_point_id'])


    initial_cond = rope.get_ivp_data()
    rope.initial_right_parts()

    if additional_info['print_system_info']:
        fixed_points_table, ropes_table, window_table, experiment_table = get_info(valid_args, rope)
        print(fixed_points_table)
        print(ropes_table)
        print(window_table)
        print(experiment_table)

    sim_instance = Simulation(rope=rope, surface=surface, surface_fill_color=surface_data['color'], 
                              rope_indexes=additional_info['rope_fragments_indexes'])
    time_data = find_time_steps(sec_num=additional_info['sec_num'], 
                            time_steps_per_second=additional_info['time_steps_per_sec'])
    sim_instance.time_step_simulation(initial_cond=initial_cond, duration=additional_info['sec_num'], 
                    time_steps_num=time_data['general_time_steps'], delay=time_data['delay'], 
                    method=additional_info['method'])
    