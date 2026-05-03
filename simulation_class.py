import numpy as np
import pygame as pg
from functions import draw_force_graph
from rope_model import RopeModel
import time



class Simulation:
    def __init__(self, rope: RopeModel, surface, surface_fill_color, rope_indexes: list = None):
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
        #for rope_indx in self.rope_indexes:
        #    self.rope.objects['rope_fragments'][rope_indx].store_force_vals = True
        print("ODE System solver running is started!")
        start_time = time.time()
        sol = self.rope.solve_ivp(initial_cond, method=method, duration=duration, fps=time_steps_num)
        end_time = time.time() - start_time
        print(f"IVP solution is found for the following time: {round(end_time, 4)} seconds!")
        
        pg.init()
        font = pg.font.SysFont("Arial", 24)

        self.rope.update_system_state(initial_cond)
        for rope_indx in self.rope_indexes:
            self.rope.objects['rope_fragments'][rope_indx].store_force_vals = True
        
        paused, running = True, True 
        t_step = 0
        
        #for t_step in range(time_steps_num):
        while running:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    running = False 
                    pg.quit()
                if event.type == pg.KEYDOWN:
                    if event.key == pg.K_SPACE:
                        paused = False if paused else True 

            updated_state = sol.y[:, t_step]
            if not paused:
                if t_step < len(sol.t) - 1:
                    t_step += 1
                else:
                    running = False
                self.rope.system(t_step, updated_state)
                self.update_coord_vals()
                self.update_speed_vals()

            self.surface.fill(self.surface_fill_color)
            self.rope.draw_model()

            if paused:
                text_surf = font.render("Press space to start the simulation!", True, (0, 0, 0))
                self.surface.blit(text_surf, (250, 550))
            pg.display.flip()     
            
            
            pg.time.delay(delay)
            

        for i in self.rope_indexes:
            force_values = self.rope.objects['rope_fragments'][i].force_vals            
            draw_force_graph(force_values, np.linspace(0, duration, len(force_values)),
                            title=f'Force (point {i})')
        draw_force_graph(self.speed_vals, time_steps=sol.t, y_label='vel(t)',
                         title=f'Velocity (weighted point)')
        draw_force_graph(self.coord_vals, time_steps=sol.t, y_label='y_coord(t)',
                         title='Vertical coordinate (weighted point)')
