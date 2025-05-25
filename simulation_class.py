import numpy as np
import matplotlib as plt
import pygame as pg
from functions import draw_force_graph


class Simulation:
    def __init__(self, rope, surface, surface_fill_color, rope_indexes: list = None):
        self.rope = rope
        self.surface = surface
        if rope_indexes:
            self.rope_indexes = rope_indexes
        else:
            self.rope_indexes = [self.rope.hor_moving_point_data['moving_points_num']]
        
        self.surface_fill_color = surface_fill_color
        self.coord_vals = []
        self.speed_vals = []
        

    def update_coord_vals(self):
        weighted_point_y_coord_val = self.rope.objects['moving_points'][-1].si_point_data['coordinates'][1]
        self.coord_vals.append(weighted_point_y_coord_val)
    
    def update_speed_vals(self):
        weighted_point_y_vel_val = self.rope.objects['moving_points'][-1].si_point_data['velocities'][1]
        self.speed_vals.append(weighted_point_y_vel_val)


    def time_step_simulation(self, initial_cond, duration, time_steps_num, delay, method='RK45'):
        for rope_indx in self.rope_indexes:
            self.rope.objects['rope_fragments'][rope_indx].store_force_vals = True
        
        sol = self.rope.solve_ivp(initial_cond, method=method, duration=duration, fps=time_steps_num)
        for t_step in range(time_steps_num):
            updated_state = sol.y[:, t_step]
            self.rope.update_system_state(updated_state)
            self.surface.fill(self.surface_fill_color)
            self.rope.draw_model()
            pg.display.flip()     
            
            self.update_coord_vals()
            self.update_speed_vals()
            pg.time.delay(delay)
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    pg.quit()

        
        for i in self.rope_indexes:
            draw_force_graph(self.rope.objects['rope_fragments'][i].force_vals)
        draw_force_graph(self.speed_vals, time_steps = sol.t, y_label = 'vel(t)')
        draw_force_graph(self.coord_vals, time_steps = sol.t, y_label = 'y_coord(t)') 
