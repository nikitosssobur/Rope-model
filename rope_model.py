import pygame as pg
import numpy as np
from points_classes import MassPoint, FixedPoint
from scipy.integrate import solve_ivp
from math import sqrt
from typing import Union


DEFAULT_REAL_METER_DIST = 100
DEFAULT_FALL_HEIGHT = 20




class RopeFragment:
    '''
    Class for rope fragment between two neighboring points!  
    '''
    def __init__(self, start_point: Union[MassPoint, FixedPoint], 
                    end_point: Union[MassPoint, FixedPoint], 
                    rest_len: Union[float, int], stiffness: Union[float, int], viscosity: Union[float, int]):
        self.rest_len = rest_len 
        self.stiffness = stiffness
        self.viscosity = viscosity
        self.start_point, self.end_point = start_point, end_point
        self.surface = None
        self.color = 'black'
        self.width = 1
        self.max_len = rest_len * 50
        self.store_force_vals = False
        self.force_vals = []


    def draw_fragment(self) -> None:
        '''
        Method for drawing rope fragment with using Pygame
        '''
        left_point_coords, right_point_coords = self.get_points_coords(mode='pix')
        pg.draw.line(self.surface, self.color, left_point_coords, right_point_coords, self.width)

    
    def get_points_coords(self, mode: str = 'si') -> tuple:
        '''
        Method returns endpoints coordinates of this RopeFragment at this moment of time
        depending on the selected mode: 
            si -- coords in SI system (meters)
            pix -- coords in pixels
        '''
        if mode == 'si':
            left_coords = self.start_point.si_point_data['coordinates']
            right_coords = self.end_point.si_point_data['coordinates']
        elif mode == 'pix': 
            left_coords = self.start_point.point_data['coordinates']
            right_coords = self.end_point.point_data['coordinates']
        return (left_coords, right_coords)
    

    def get_points_vels(self) -> tuple:
        '''
        Method returns endpoints velocities in SI system (meters per second)
        '''
        left_vels = self.start_point.si_point_data['velocities']
        right_vels = self.end_point.si_point_data['velocities']
        return (left_vels, right_vels)


    def vector_coords(self, opposite_vector: bool = False) -> np.ndarray:
        '''
        Coordinates of the RopeFragment vector
        '''
        left_coords, right_coords = self.get_points_coords()
        if opposite_vector:
            vector_coords = [left_coords[0] - right_coords[0], left_coords[1] - right_coords[1]]
        else:
            vector_coords = [right_coords[0] - left_coords[0], right_coords[1] - left_coords[1]]
        return np.array(vector_coords)


    def vel_vector_coords(self, opposite_vector: bool = False) -> np.ndarray:
        '''
        Velocities of the RopeFragment vector
        '''
        left_vels, right_vels = self.get_points_vels()
        if opposite_vector:
            vector_vels = [left_vels[0] - right_vels[0], left_vels[1] - right_vels[1]]
        else:
            vector_vels = [right_vels[0] - left_vels[0], right_vels[1] - left_vels[1]]
        return np.array(vector_vels)


    def rope_fragment_len(self):
        '''
        RopeFragment lenght at this moment in time
        '''
        vector_coords = self.vector_coords()
        dist = np.linalg.norm(vector_coords)
        if dist >= self.max_len:
            raise ValueError("Rope lenght can't be greater than maximum length!")
        return dist


    def force(self, main_point: str = 'left') -> np.ndarray:
        '''
        Tension force in the RopeFragment at this moment of time derived using 
        Hooke's law and Newton's law of viscous friction (this force is used 
        as a component of the right part of the system of ODEs in RopeModel class)
        '''
        if main_point == 'left':
            vector, vel_vector = self.vector_coords(), self.vel_vector_coords()
        elif main_point == 'right':
            vector, vel_vector = self.vector_coords(opposite_vector = True), self.vel_vector_coords(opposite_vector = True)
        else:
            raise ValueError("Incorrect value of main_point parameter! Must be 'left' or 'right'.")

        
        vector_norm = self.rope_fragment_len()
        if vector_norm <= self.rest_len:
            self.force_val = 0
            if self.store_force_vals:
                self.force_vals.append(self.force_val)
            return np.array([0.0, 0.0])
        
        
        vector /= vector_norm
        elastic = self.stiffness * (vector_norm - self.rest_len)
        viscous = self.viscosity * np.dot(vector, vel_vector)
        self.force_vec = (elastic + viscous) * vector / self.rest_len
        self.force_val = np.linalg.norm(self.force_vec)
        if self.store_force_vals:
            self.force_vals.append(self.force_val)
        return self.force_vec
    


class RopeModel:
    '''
    General rope model (contains method for horizontal rope and vertical rope)
    '''
    
    def __init__(self, surface, hor_rope_data, real_meter_dist=None, fall_height=None):
        self.surface = surface
        self.fixed_point_data = hor_rope_data['fixed_point_data']
        self.hor_moving_point_data = hor_rope_data['moving_point_data'] 
        self.hor_rope_mass = hor_rope_data['rope_mass']
        self.hor_rope_stiffness = hor_rope_data['rope_stiffness']
        self.hor_rope_viscosity = hor_rope_data['rope_viscosity']
        self.hor_rope_tension = hor_rope_data['tension_force']
        self.fall_height = fall_height if fall_height is not None and fall_height >= 0 else DEFAULT_FALL_HEIGHT
        self.real_meter_dist = real_meter_dist if real_meter_dist is not None and real_meter_dist >= 0 else DEFAULT_REAL_METER_DIST
        self.set_convert_coeff(real_meter_dist=self.real_meter_dist)
        self.set_hor_rest_len()
        self.objects = {'fixed_points': [], 
                        'rope_fragments':[],
                        'moving_points':[],
                        'surface': self.surface}
        

    def create_fixed_points(self):
        for i in range(len(self.fixed_point_data['x_coords'])):
            fixed_point = FixedPoint(self.fixed_point_data['x_coords'][i], 
                                self.fixed_point_data['y_coords'][i], 
                                self.fixed_point_data['radius'])
            disp_vec = fixed_point.point_data['coordinates'] if i == 0 else disp_vec
            fixed_point.set_convert_data(scale = self.pix_per_metr, displace_vector = disp_vec)
            fixed_point.set_surface(self.surface)
            self.objects['fixed_points'].append(fixed_point)

    
    def create_free_end_rope(self, attach_point_id: int, points_velocities = None):
        self.attach_point_id = attach_point_id
        top_point = self.objects['moving_points'][self.attach_point_id]
        points_num = self.ver_rope_point_num
        x_coord = top_point.point_data['coordinates'][0]
        fragment_rest_len = self.get_rope_fragment_rest_len(rope_type = 'ver')
        start_coord = top_point.point_data['coordinates'][1]
        end_coord = start_coord + self.ver_rope_rest_len * self.pix_per_metr
        vertical_coords = np.linspace(start_coord, end_coord, points_num + 1)[1:]
        vertical_coords = [(x_coord, y_coord) for y_coord in vertical_coords]
        stiffness = self.ver_rope_stiffness        
        viscosity = self.ver_rope_viscosity        
        mass = self.get_moving_point_mass(rope_type = 'ver')
        radius = self.ver_rope_point_radius
        if points_velocities is None:
            y_vels = self.get_linear_vels_distribution(points_num + 1, 0, sqrt(2 * 9.81 * self.fall_height))[1:]
            points_velocities = [(0, y_vel * self.pix_per_metr) for y_vel in y_vels]

        for i in range(points_num):
            if i == points_num - 1:
                mass = self.weighted_point_mass
                radius = self.weighted_point_radius

            if i != 0:
                top_point = self.objects['moving_points'][-1] 
            
            self.create_moving_point(id = len(self.objects['moving_points']), 
                                    mass = mass, radius = radius, 
                                    coords = vertical_coords[i], 
                                    velocities = points_velocities[i])
            bottom_point = self.objects['moving_points'][-1]
            self.create_rope_fragment(top_point, bottom_point, stiffness, viscosity, fragment_rest_len)


        self.general_points_num = len(self.objects['moving_points'])
    

    def create_rope_fragment(self, left_endpoint, right_endpoint, stiffness = None, viscosity = None, rest_len = None):
        if rest_len is None:
            rest_len = self.get_rope_fragment_rest_len()

        if stiffness is None:
            stiffness = self.hor_rope_stiffness       
                                                        

        if viscosity is None:
            viscosity = self.hor_rope_viscosity          
                                                        

        rope_fragment = RopeFragment(left_endpoint, right_endpoint, rest_len, stiffness, viscosity)
        rope_fragment.surface = self.surface
        self.objects['rope_fragments'].append(rope_fragment)

        
    def get_linear_vels_distribution(self, point_num: int, start_vel: float, end_vel: float) -> np.ndarray:
        vels = np.linspace(start_vel, end_vel, point_num)
        return vels


    def create_moving_point(self, id: int, mass: float, radius: float, coords, velocities) -> None:
        moving_point = MassPoint(id, mass, radius, (coords[0], coords[1]), 
                                (velocities[0], velocities[1]))
        moving_point.set_surface(self.surface)
        moving_point.set_convert_data(scale = self.pix_per_metr, 
                                      displace_vector = self.objects['fixed_points'][0].point_data['coordinates'])
        self.objects['moving_points'].append(moving_point)
        

    def create_model(self):
        self.create_fixed_points()
        hor_coords = np.linspace(self.objects['fixed_points'][0].x, 
                                 self.objects['fixed_points'][1].x, 
                                 self.hor_moving_point_data['moving_points_num'] + 2)[1 :-1 :]
        
        vertical_coord = self.fixed_point_data['y_coords'][0]
        self.set_moving_point_mass(self.hor_rope_mass, self.hor_moving_point_data['moving_points_num'])
        
        for i in range(self.hor_moving_point_data['moving_points_num']):
            self.create_moving_point(i, self.hor_moving_point_mass, 
                                     self.hor_moving_point_data['moving_points_radius'],
                                     (hor_coords[i], vertical_coord), 
                                     self.hor_moving_point_data['moving_points_velocities'])
            

        left_point = self.objects['fixed_points'][0] 
        right_point = self.objects['moving_points'][0]
        self.create_rope_fragment(left_point, right_point)
        
        for j in range(self.hor_moving_point_data['moving_points_num'] - 1):
            left_point = self.objects['moving_points'][j]
            right_point = self.objects['moving_points'][j + 1]
            self.create_rope_fragment(left_point, right_point) 

        left_point = self.objects['moving_points'][-1] 
        right_point = self.objects['fixed_points'][1]
        self.create_rope_fragment(left_point, right_point)              

        
    def draw_model(self):
        for fixed_point in self.objects['fixed_points']:
            fixed_point.draw_point()
        

        for moving_point in self.objects['moving_points']:
            moving_point.draw_point()


        for rope_fragment in self.objects['rope_fragments']:
            rope_fragment.draw_fragment()


    def set_moving_point_mass(self, rope_mass: float, moving_points_num: int, rope_type = 'hor'):
        if rope_type == 'hor':
            self.hor_moving_point_mass = rope_mass / (moving_points_num) 
        elif rope_type == 'ver':
            self.ver_moving_point_mass = rope_mass / (moving_points_num - 1) 
        
    
    def get_moving_point_mass(self, rope_type = 'hor'):
        if rope_type == 'hor': return self.hor_moving_point_mass
        elif rope_type == 'ver': return self.ver_moving_point_mass


    def add_mass(self, mass: float, point_indx: int):
        masspoint = self.objects['moving_points'][point_indx]
        masspoint.mass += mass if mass >= 0 else 0
    

    def get_fixed_points_distanse(self):
        x_coords = self.fixed_point_data['x_coords']
        dist = x_coords[1] - x_coords[0]
        return dist


    def set_hor_rest_len(self):
        stretch_rope_len = self.get_fixed_points_distanse() / self.pix_per_metr
        hor_rest_len = stretch_rope_len / (self.hor_rope_tension / self.hor_rope_stiffness + 1)
        if hor_rest_len > 0:
            self.hor_rest_len = hor_rest_len
        else:
            print(f"""Horizontal rest len: {hor_rest_len},\n 
                    Stretched  rope len: {stretch_rope_len},\n 
                    Horizontal rope stiffness: {self.hor_rope_stiffness}""")
            raise ValueError("Rest len of the horizontal rope must be a positive float!")


    def set_convert_coeff(self, real_meter_dist = DEFAULT_REAL_METER_DIST):
        dist = self.get_fixed_points_distanse()
        self.pix_per_metr = dist / real_meter_dist


    def get_rope_fragment_rest_len(self, rope_type = 'hor'):
        if rope_type == 'hor':
            rope_fragments_num = self.hor_moving_point_data['moving_points_num'] + 1
            hor_fragment_rest_len = self.hor_rest_len / rope_fragments_num
            return hor_fragment_rest_len
        
        elif rope_type == 'ver':
            rope_fragments_num = self.ver_rope_point_num 
            ver_fragment_rest_len = self.ver_rope_rest_len / rope_fragments_num
            return ver_fragment_rest_len
        

    def set_ver_rope_data(self, ver_rope_data):
        self.ver_rope_rest_len = ver_rope_data['rest_len'] / self.pix_per_metr
        self.ver_rope_mass = ver_rope_data['rope_mass']
        self.ver_rope_stiffness = ver_rope_data['rope_stiffness']
        self.ver_rope_viscosity = ver_rope_data['rope_viscosity']
        self.ver_rope_point_num = ver_rope_data['point_num'] 
        self.ver_rope_point_radius = ver_rope_data['ver_rope_point_radius']
        self.weighted_point_mass = ver_rope_data['weighted_point_mass'] 
        self.weighted_point_radius = ver_rope_data['weighted_point_radius'] 
        self.set_moving_point_mass(self.ver_rope_mass, self.ver_rope_point_num, rope_type = 'ver')


    def split_data(self, data_list):
        '''
        Splitting the data_list containing info about coordinates and velocities of each point in one 
        general 1d-list into the following format:
        [data_1, data_2, ..., data_point_num], where data_i is a data list [x, dx, y, dy] 
        for i-th point of the system.
        '''
        splitted_data = []
        for i in range(0, len(data_list), 4):
            splitted_data.append(data_list[i : i + 4])

        return splitted_data

    
    def system(self, t, data_list):
        '''
        A method for solving the system of differential equations. 
        This function is passed as an argument into solve_ivp from Scipy. 
        '''
        if self.first_num_step:
            self.first_num_step = False
            return self.init_right_parts
        

        self.update_system_state(data_list)
        splitted_data = self.split_data(data_list)
        updated_right_parts = []  
        for point_num in range(self.general_points_num):
            '''
            This loop generates 4 diff equations for each point of the system:
            2 equations for x and x' (x coordinate and x velocity (the first derivative of x by time))
            2 equations for y and y'
            '''
            point_mass = self.objects['moving_points'][point_num].point_data['mass'] 
            if point_num <= self.hor_moving_point_data['moving_points_num'] - 1:
                #Horizontal rope points
                left_rope = self.objects['rope_fragments'][point_num]
                right_rope = self.objects['rope_fragments'][point_num + 1]
                #right_rope.store_force_vals = True if np.any(np.isclose(t, self.t_eval, atol=2.5*1e-3, rtol=1e-5)) else False   #atol=1.747*1e-4  #atol=5.4*1e-5, rtol=1e-5                     
                

                if point_num == self.attach_point_id:
                    #Attachment point case
                    indx = self.hor_moving_point_data['moving_points_num'] + 1  
                    vert_rope = self.objects['rope_fragments'][indx]


            elif self.hor_moving_point_data['moving_points_num'] - 1 < point_num:
                #Vertical rope points case 
                top_rope = self.objects['rope_fragments'][point_num + 1]
                if point_num < self.general_points_num - 1:
                    #Vertical rope points (case of unweighted point)
                    down_rope = self.objects['rope_fragments'][point_num + 2]
            
            '''
            Definitions:
            w  = x -- x coord of the point in the previous time step
            dw = w' = z = x' -- velocity of the x coord
            dz = z' = x'' -- acceleration of the x coord  
            u  = y -- y coord of the point in the previous time step
            du = u' = v = y' -- velocity of the y coord 
            dv = v' = y'' -- acceleration of the y coord 
            '''
            w, z, u, v = splitted_data[point_num]
            
            if point_num < self.general_points_num - 1: 
                if point_num == self.attach_point_id: #attachment point case
                    left_rope_force = left_rope.force(main_point = 'right')
                    right_rope_force = right_rope.force()
                    if self.general_points_num == self.hor_moving_point_data['moving_points_num']:
                        vert_rope_force = np.array([0.0, 0.0])
                    else:
                        vert_rope_force = vert_rope.force()
                    dz = (left_rope_force[0] + right_rope_force[0] + vert_rope_force[0]) / point_mass    #+ 9.81
                    dv = (left_rope_force[1] + right_rope_force[1] + vert_rope_force[1]) / point_mass     #+ 9.81 


                elif point_num > self.hor_moving_point_data['moving_points_num'] - 1:
                    top_rope_force = top_rope.force(main_point = 'right')
                    down_rope_force = down_rope.force()
                    dz = (top_rope_force[0] + down_rope_force[0]) / point_mass    #+ 9.81
                    dv = (top_rope_force[1] + down_rope_force[1]) / point_mass + 9.81 


                else:
                    left_rope_force = left_rope.force(main_point = 'right')
                    right_rope_force = right_rope.force()
                    dz = (left_rope_force[0] + right_rope_force[0]) / point_mass   #+ 9.81
                    dv = (left_rope_force[1] + right_rope_force[1]) / point_mass #+ 9.81 
               

            else: #weighted point case
                
                top_rope_force = top_rope.force(main_point = 'right')
                dz = top_rope_force[0] / point_mass #+ 9.81
                dv = top_rope_force[1] / point_mass + 9.81
            #Diff eqs for mass point with number point_num

            dw = z
            du = v
            updated_right_parts += [dw, dz, du, dv]

        return updated_right_parts


    def get_ivp_data(self):
        '''
        Inverse function (in some sence) to split_data function. 
        This function converts coordinates and velocities of each separate point into the general 
        list, which contains all coords and vels.
        '''
        data = []
        for point in self.objects['moving_points']:
            coords = point.si_point_data['coordinates']
            vels = point.si_point_data['velocities']
            data += [coords[0], vels[0], coords[1], vels[1]]
        return data
    

    def initial_right_parts(self):
        self.init_right_parts = []
        self.first_num_step = True
        for i, point in enumerate(self.objects['moving_points']):
            mass, vels = point.point_data['mass'], point.si_point_data['velocities']
            if i < self.hor_moving_point_data['moving_points_num']: 
                right_part_point_data = [vels[0], self.hor_rope_tension / mass, vels[1], 0]
            else:
                right_part_point_data = [vels[0], 0, vels[1], 9.81]
            self.init_right_parts += right_part_point_data


    def update_system_state(self, new_state_list):
        state = self.split_data(new_state_list)
        for i, point in enumerate(self.objects['moving_points']):
            new_point_state = state[i]
            point.update_coords(new_point_state[0], new_point_state[2])  #/ self.scale, new_point_state[2] / self.scale)
            point.update_velocities(new_point_state[1], new_point_state[3])


    def solve_ivp(self, initial_cond, duration = 1, start_time = 0, fps = 60, method = 'RK45'): #RK45 
        '''
        Method for solving a system of differential equations using 
        a numerical method at each step of the game cycle  
        '''
        fps = fps
        duration = duration # 1 sec
        t_span = (start_time, start_time + duration) #(0, 1)
        self.t_eval = np.linspace(t_span[0], t_span[1], fps) #(0, 0.1, 0.2, 1) 
        solution = solve_ivp(self.system, method = method, t_span = t_span, y0 = initial_cond, 
                             t_eval = self.t_eval)
        return solution


    def set_weighted_point_data(self, radius, mass):
        weighted_point = self.objects['moving_points'][-1]
        weighted_point.set_mass(mass)
        weighted_point.set_radius(radius)