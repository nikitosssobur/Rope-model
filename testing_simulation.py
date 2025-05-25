import pygame as pg
from simulation_class import Simulation
from rope_model import RopeModel
from functions import unpack_model_config, find_time_steps


hor_rope_data, ver_rope_data, surface_data, additional_info = unpack_model_config(r"./default_model_config.json")

pg.init()
surface = pg.display.set_mode((surface_data['width'],  surface_data['height']))

rope = RopeModel(surface, hor_rope_data=hor_rope_data)
rope.set_ver_rope_data(ver_rope_data)
rope.create_model()
rope.create_free_end_rope(attach_point_id=additional_info['attach_point_id'])


initial_cond = rope.get_ivp_data()
rope.initial_right_parts()

sim_instance = Simulation(rope=rope, surface=surface, surface_fill_color=surface_data['color'])
time_data = find_time_steps(sec_num=additional_info['sec_num'], 
                            time_steps_per_second=additional_info['time_steps_per_sec'])
sim_instance.time_step_simulation(initial_cond=initial_cond, duration=additional_info['sec_num'], 
                    time_steps_num=time_data['general_time_steps'], delay=time_data['delay'], 
                    method='RK23')#additional_info['method'])



