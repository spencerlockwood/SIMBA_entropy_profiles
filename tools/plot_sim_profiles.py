## Simulation plotting functions
## Load libraries

import warnings
warnings.filterwarnings('ignore')

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colorbar import make_axes, make_axes_gridspec

import numpy as np
import pynbody as pnb

import periodictable as pt

import copy

import os

import astropy.constants as c
import astropy.units as u
import astropy.io.fits as pyfits

from num2tex import num2tex, configure
configure(exp_format='times', display_singleton=True)

import logging

from palettable.cartocolors.sequential import DarkMint_7, BurgYl_7, Purp_7, Sunset_7, agGrnYl_7


from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units, G#, median_R_values, #Rxxx_values, Rxxx_conversions
import group_xray_profiles_observations_plotting_functions as xobs
# from oppenheimer21.groups_profiles_all import run_plotting_script as opp21_plot







# class CustomString:
#     def __init__(self, value):
#         self.value = value

#     def __str__(self):
#         return self.value

#     def __repr__(self):
#         return f"CustomString({self.value!r})"

#     def __eq__(self, other):
#         # Ensure that equality is only based on object identity, not content
#         return self is other

#     def __hash__(self):
#         # Use the object's id for a unique hash
#         return id(self)
    
#     def __len__(self):
#         return len(self.value)






def set_list_duplicates_to_none(lst, sort_list=None, 
                                all_duplicates=True, specific_value=None, 
                                remove_followers=True, return_sort_list=False):
    lst_copy = lst.copy()
    if sort_list is None:
        sort_lst = lst.copy()
    else:
        sort_lst = sort_list.copy()
    
    idx = 0
    incr = 1
    stopidx = len(lst_copy)
    if remove_followers:
        idx = stopidx - 1
        stopidx = incr = -1
        
    allow_to_replace = True
        
    while idx != stopidx:
        elem = sort_lst[idx]
    
        if not all_duplicates and elem!=specific_value:
            allow_to_replace = False
    
        if sort_lst.count(elem) > 1 and allow_to_replace:
            lst_copy[idx] = None
            sort_lst[idx] = None
            
        idx += incr
        allow_to_replace = True
        
    if return_sort_list:
        return lst_copy, sort_lst
    else:
        return lst_copy

    
    
    
    

# def add_cax_to_figure(fig, nrows=1, ncols=1, location='right'):
#     """Add a colorbar to the figure at the specified location after subplots are created."""
#     # Define a new GridSpec for the figure with an extra row or column for the colorbar
#     if location == 'right':
#         gs = gridspec.GridSpec(nrows, ncols+1, figure=fig)  # Extra column for colorbar
#         cax = fig.add_subplot(gs[:, -1])  # Colorbar in the last column
#     elif location == 'left':
#         gs = gridspec.GridSpec(nrows, ncols+1, figure=fig)  # Extra column for colorbar
#         cax = fig.add_subplot(gs[:, 0])  # Colorbar in the first column
#     elif location == 'top':
#         gs = gridspec.GridSpec(nrows+1, ncols, figure=fig)  # Extra row for colorbar
#         cax = fig.add_subplot(gs[0, :])  # Colorbar in the first row
#     elif location == 'bottom':
#         gs = gridspec.GridSpec(nrows+1, ncols, figure=fig)  # Extra row for colorbar
#         cax = fig.add_subplot(gs[-1, :])  # Colorbar in the last row
#     else:
#         raise ValueError("Invalid location. Choose from 'left', 'right', 'top', 'bottom'.")
        
#     cax.set_visible(False)
        
#     return cax

#     # Create the colorbar
# #     return plt.colorbar(img, cax=cax, orientation='horizontal' if location in ['top', 'bottom'] else 'vertical')
    
    
    
def add_cax_to_figure(fig, ax, ax_is_fig=False, x=0, y=0, width=0, height=0, 
                      x_in_inches=False, y_in_inches=False, width_in_inches=False, height_in_inches=False):
    """Add a colorbar to the figure at the specified location after subplots are created."""
    
    width_inches, height_inches = fig.get_size_inches()
    if x_in_inches:
        x /= width_inches
    if width_in_inches:
        width /= width_inches    
    if y_in_inches:
        y /= height_inches
    if height_in_inches:
        height /= height_inches
        
#     else:
#         width_inches, height_inches = 1, 1
        
#     print()
#     print('width_inches:', width_inches)
#     print('height_inches:', height_inches)
#     print()
    
    if ax_is_fig:
        x0, x1, y0, y1 = 0, 1, 0, 1
    else:
        ax_pos = ax.get_position()
        x0, x1, y0, y1 = ax_pos.x0, ax_pos.x1, ax_pos.y0, ax_pos.y1
    
#     x0 *= width_inches
#     x1 *= width_inches
#     y0 *= height_inches
#     y1 *= height_inches
    
    ax_width = x1 - x0
    ax_height = y1 - y0
    
    x_norm = x0 + x*ax_width
    y_norm = y0 + y*ax_height
    
    width_norm = width*ax_width
    height_norm = height*ax_height
    
    left_norm = x_norm - 0.5*width_norm
    bottom_norm = y_norm - 0.5*height_norm
    
#     if in_inches:
#         left_norm /= width_inches
#         bottom_norm /= height_norm
#         width_norm /= width_inches
#         height_norm /= height_inches
    
    cbar_pos = [left_norm, bottom_norm, width_norm, height_norm]
    cax = fig.add_axes(cbar_pos, visible=True)
        
    return cax
    
    
    


## Figure out correct errors to use on individual profiles
def calc_poisson_err(x, N, log=False):
    if log:
#         return np.log10((10**x)/(N**0.5))
        return x/N
    else:
        return x/(N**0.5)


# def plot_extra_lines(axis=plt, xunits=None, plot_kwargs=None):
    
#     if plot_kwargs['solar_line']['plot']:
#         if plot_kwargs['solar_line']['show_labels']:
#             labels = {
# #                 'Zsol':r'$Z_{\odot}$',
#                 'Zsol':'#lines##Solar',
#             }
#         else:
#             labels = {
#                 'Zsol':None,
#             }
            
#         if not plot_kwargs['yaxis']['is_log']:
#             axis.axhline(1, color=plot_kwargs['solar_line']['color'], 
#                          zorder=plot_kwargs['solar_line']['zorder'], 
#                          ls=plot_kwargs['solar_line']['ls'], lw=plot_kwargs['solar_line']['lw'], 
#                          label=labels['Zsol'])
            
#         elif plot_kwargs['yaxis']['is_log']:
#             axis.axhline(np.log10(1), color=plot_kwargs['solar_line']['color'], 
#                          zorder=plot_kwargs['solar_line']['zorder'], 
#                          ls=plot_kwargs['solar_line']['ls'], lw=plot_kwargs['solar_line']['lw'], 
#                          label=labels['Zsol'])
    
    
#     if plot_kwargs['Rxxx_lines']['plot']:
#         if plot_kwargs['Rxxx_lines']['show_labels']:
#             labels = {
#                 'R2500':r'#lines##$R_{2500}$',
#                 'R500':r'#lines##$R_{500}$',
#                 'R200':r'#lines##$R_{200}$',
#             }
#         else:
#             labels = {
#                 'R2500':None,
#                 'R500':None,
#                 'R200':None,
#             }
            
# #         if 'physical' not in plot_kwargs['xaxis']['units'].lower():
# #             main_line = 1 if not plot_kwargs['xaxis']['is_log'] else np.log10(1)
# #         else:
# #             main_line = Rxxx_values[ if not plot_kwargs['xaxis']['is_log'] else np.log10(1)
        
#         if not plot_kwargs['xaxis']['is_log']:
# #             axis.axvline(1, color=plot_kwargs['Rxxx_lines']['color'], 
# #                          zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
# #                          label=labels[plot_kwargs['xaxis']['units'].upper()])

#             if 'r200' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(1, color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels[plot_kwargs['xaxis']['units'].upper()])
#                 axis.axvline(Rxxx_conversions['R500/R200'], color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R500'])
#                 axis.axvline(Rxxx_conversions['R2500/R200'], color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-.', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R2500'])

#             elif 'r500' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(1, color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels[plot_kwargs['xaxis']['units'].upper()])
#                 axis.axvline(Rxxx_conversions['R200/R500'], color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R200'])
#                 axis.axvline(Rxxx_conversions['R2500/R500'], color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-.', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R2500'])

#             elif 'r2500' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(1, color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels[plot_kwargs['xaxis']['units'].upper()])
#                 axis.axvline(Rxxx_conversions['R200/R2500'], color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R200'])
#                 axis.axvline(Rxxx_conversions['R500/R2500'], color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-.', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R500'])
                
#             elif 'physical' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(Rxxx_values['R500'].in_units(xunits), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R500'])
#                 axis.axvline(Rxxx_values['R200'].in_units(xunits), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R200'])
#                 axis.axvline(Rxxx_values['R2500'].in_units(xunits), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R2500'])
                

#         elif plot_kwargs['xaxis']['is_log']:
# #             axis.axvline(np.log10(1), color=plot_kwargs['Rxxx_lines']['color'], 
# #                          zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
# #                          label=labels[plot_kwargs['xaxis']['units'].upper()])

#             if 'r200' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(np.log10(1), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels[plot_kwargs['xaxis']['units'].upper()])
#                 axis.axvline(np.log10(Rxxx_conversions['R500/R200']), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R500'])
#                 axis.axvline(np.log10(Rxxx_conversions['R2500/R200']), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-.', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R2500'])

#             elif 'r500' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(np.log10(1), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels[plot_kwargs['xaxis']['units'].upper()])
#                 axis.axvline(np.log10(Rxxx_conversions['R200/R500']), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R200'])
#                 axis.axvline(np.log10(Rxxx_conversions['R2500/R500']), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-.', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R2500'])

#             elif 'r2500' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(np.log10(1), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels[plot_kwargs['xaxis']['units'].upper()])
#                 axis.axvline(np.log10(Rxxx_conversions['R200/R2500']), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R200'])
#                 axis.axvline(np.log10(Rxxx_conversions['R500/R2500']), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-.', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R500'])
                
#             elif 'physical' in plot_kwargs['xaxis']['units'].lower():
#                 axis.axvline(np.log10(Rxxx_values['R500'].in_units(xunits)), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='-', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R500'])
#                 axis.axvline(np.log10(Rxxx_values['R200'].in_units(xunits)), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R200'])
#                 axis.axvline(np.log10(Rxxx_values['R2500'].in_units(xunits)), color=plot_kwargs['Rxxx_lines']['color'], 
#                              zorder=plot_kwargs['Rxxx_lines']['zorder'], ls='--', lw=plot_kwargs['Rxxx_lines']['lw'], 
#                              label=labels['R2500'])
                
                
                
def plot_extra_lines_new(labels, orientations, values, ls, lw, colors, zorders, axis=plt):# plot_kwargs=None):
#     for line_idx in range(len(plot_kwargs['labels'])):
#         if plot_kwargs['orientation'][line_idx].lower() == 'v':
#             axis.axvline(plot_kwargs['values'][line_idx], 
#                          ls=plot_kwargs['ls'][line_idx], 
#                          lw=plot_kwargs['lw'][line_idx], 
#                          color=plot_kwargs['color'][line_idx], 
#                          zorder=plot_kwargs['zorder'][line_idx], 
#                          label=plot_kwargs['labels'][line_idx])
#         elif plot_kwargs['orientation'][line_idx].lower() == 'h':
#             axis.axhline(plot_kwargs['values'][line_idx], 
#                          ls=plot_kwargs['ls'][line_idx], 
#                          lw=plot_kwargs['lw'][line_idx], 
#                          color=plot_kwargs['color'][line_idx], 
#                          zorder=plot_kwargs['zorder'][line_idx], 
#                          label=plot_kwargs['labels'][line_idx])
#         else:
#             raise Exception('Orientation of extra lines must be either vertical (v) or horizontal (h).')
            
            
    for line_idx in range(len(values)):
        if orientations[line_idx].lower() == 'v':
            axis.axvline(values[line_idx], 
                         ls=ls[line_idx], 
                         lw=lw[line_idx], 
                         color=colors[line_idx], 
                         zorder=zorders[line_idx], 
                         label=labels[line_idx])
        elif orientations[line_idx].lower() == 'h':
            axis.axhline(values[line_idx], 
                         ls=ls[line_idx], 
                         lw=lw[line_idx], 
                         color=colors[line_idx], 
                         zorder=zorders[line_idx], 
                         label=labels[line_idx])
        else:
            raise Exception('Orientation of extra lines must be either vertical (v) or horizontal (h).')
                
                

                
def plot_individual_profiles(profiles_info, code, prop, bin_, type_='regular', classification_='all',
                             axis=plt, ind_color='black', ind_ls='-', sm=None,
                             metal_idx=None, metal_idx_num=None, metal_idx_denom=None,
                             og_prop_yunits=None, prop_yunits=None, label=None,
                             plot_kwargs=None, fig_kwargs=None, 
                             current_axis_info=None):
    
    ind_plot_color = ind_color
        
    res_number_counts_ = {}
    for halo in profiles_info[code]['profiles'][bin_][classification_+'_id']:
        if type_.lower() == 'regular':
            bin_prop_val_ = profiles_info[code]['profiles'][bin_][classification_+'_profiles'][halo]['y'][prop]
        elif type_.lower() == 'metallicity':
            bin_prop_val_ = profiles_info[code]['profiles'][bin_][classification_+'_profiles'][halo]['y'][prop][:,metal_idx]
        elif type_.lower() == 'metallicity ratio':
            bin_prop_val_ = profiles_info[code]['profiles'][bin_][classification_+'_profiles'][halo]['y'][prop][:,metal_idx_num,metal_idx_denom]
        else:
            raise Exception('type_ must be regular, metallicity, or metallicity ratio')

        ## Switch to correct units
        try:
            bin_prop_val_ = bin_prop_val_.in_units(prop_yunits)
        except:
            bin_prop_val_ = pnb.array.SimArray(bin_prop_val_, units=og_prop_yunits)
            bin_prop_val_ = bin_prop_val_.in_units(prop_yunits)
        

        if plot_kwargs['plot_ind']:
            if fig_kwargs['cbar']['plot']:
                halo_idx = list(profiles_info[code]['profiles'][bin_][classification_+'_id']).index(halo)
                bin_val_ = profiles_info[code]['profiles'][bin_]['halo_info'][fig_kwargs['binning']['prop'][0]][fig_kwargs['binning']['prop'][1]][halo_idx]
                if fig_kwargs['cbar']['log_cbar']:
                    bin_val_ = np.log10(bin_val_)
                ind_cmap_color_ = fig_kwargs['cbar']['cm']((bin_val_-sm.get_clim()[0])/(sm.get_clim()[1]-sm.get_clim()[0]))
                ind_plot_color = ind_cmap_color_

            axis.plot(profiles_info[code]['XAXIS'], bin_prop_val_, 
                      color=ind_plot_color, ls=ind_ls, lw=plot_kwargs['lw'], 
                      alpha=plot_kwargs['alpha'], drawstyle=plot_kwargs['drawstyle'], 
                      zorder=plot_kwargs['zorder'])
            if plot_kwargs['plot_poisson_err']:
                poisson_err_ = calc_poisson_err(bin_prop_val_, 
                                                profiles_info[code]['profiles'][bin_][classification_+'_profiles'][halo]['y']['n'], 
                                                log=fig_kwargs['yaxis']['is_log'])
                axis.fill_between(profiles_info[code]['XAXIS'], 
                                      bin_prop_val_ + poisson_err_, bin_prop_val_ - poisson_err_, 
                                      color=ind_plot_color, alpha=plot_kwargs['alpha'], 
                                      zorder=plot_kwargs['zorder'])

        if plot_kwargs['calc_number_lines']:
            x_res_ = plot_resolution_lines(profiles_info, code, bin_, halo, classification_=classification_, 
                                           axis=axis, color=ind_plot_color, plot_kwargs=plot_kwargs)

        if plot_kwargs['calc_number_lines']:
            if x_res_ in res_number_counts_.keys():
                res_number_counts_[x_res_] += 1
            else:
                res_number_counts_[x_res_] = 1
    

    if plot_kwargs['calc_number_lines']:
        all_x_res_values = []
        for x_, count_ in res_number_counts_.items():
            if plot_kwargs['plot_ind_number_lines']:
                axis.text(s=count_, x=1.01*x_, y=0.99*max(axis.get_ylim()), va='top', ha='left', color=ind_plot_color)
            if plot_kwargs['plot_avg_number_lines']:
                all_x_res_values += [float(x_)]*count_
                
        if plot_kwargs['plot_avg_number_lines']:
            all_x_res_values = np.array(all_x_res_values)
#             median_x_res_value = np.nanmedian(all_x_res_values)
            mean_x_res_value = np.nanmean(all_x_res_values)
#             std_x_res_value = np.nanstd(all_x_res_values)
            percentile_x_res_values = np.nanquantile(all_x_res_values, [0.16, 0.84])
            axis.axvline(x=mean_x_res_value, color=ind_color, 
                         ls=ind_ls, alpha=plot_kwargs['number_line_alpha'], 
                         lw=plot_kwargs['number_line_lw'], zorder=plot_kwargs['number_line_zorder']+1)#, 
#                          label=label+' (ind)')
            axis.fill_between(x=[percentile_x_res_values[0], percentile_x_res_values[1]], 
                              y1=fig_kwargs['yaxis']['lo_lims'][current_axis_info['row']] if fig_kwargs['yaxis']['set_lims'] 
                              else min(axis.get_ylim()), 
                              y2=fig_kwargs['yaxis']['hi_lims'][current_axis_info['row']] if fig_kwargs['yaxis']['set_lims'] 
                              else max(axis.get_ylim()), 
                              color=ind_color, alpha=0.2, zorder=plot_kwargs['number_line_zorder']+1)
#             axis.fill_between(x=[mean_x_res_value-std_x_res_value, mean_x_res_value+std_x_res_value], 
#                               y1=fig_kwargs['yaxis']['lo_lims'][current_axis_info['row']] if fig_kwargs['yaxis']['set_lims'] 
#                               else min(axis.get_ylim()), 
#                               y2=fig_kwargs['yaxis']['hi_lims'][current_axis_info['row']] if fig_kwargs['yaxis']['set_lims'] 
#                               else max(axis.get_ylim()), 
#                               color=ind_color, alpha=0.2, zorder=plot_kwargs['number_line_zorder']+1)
            


def plot_resolution_lines(profiles_info, code, bin_, halo, classification_='all', axis=plt, color='black', plot_kwargs=None):
    ## Finds all radial bins that have >= N number of particles
    total_counts = profiles_info[code]['profiles'][bin_][classification_+'_profiles'][halo]['y']['n']
    resolved_counts = total_counts[total_counts >= plot_kwargs['res_limit']]
    resolved_idxs = np.where(total_counts >= plot_kwargs['res_limit'])[0]
    unresolved_idxs = np.where(total_counts < plot_kwargs['res_limit'])[0]
    
    if plot_kwargs['res_limit_method'].lower() == 'first_bin' and len(unresolved_idxs) > 0:
        binnum_resolution = int(np.max(unresolved_idxs)) ## first radial bin (going inwards) that has < N particles
    elif plot_kwargs['res_limit_method'].lower() == 'last_bin' and len(resolved_idxs) > 0:
        binnum_resolution = int(np.min(resolved_idxs)) ## last radial bin (going inwards) that has >= N particles
    elif plot_kwargs['res_limit_method'].lower() == 'min_particles' and len(resolved_idxs) > 0:
        try:
            min_idx = np.argmin(resolved_counts)[0]
        except:
            min_idx = np.argmin(resolved_counts)
        binnum_resolution = int(resolved_idxs[min_idx]) ## Radial bin that has >= N and min # of particles
    else:
        binnum_resolution = 0

    x_resolution = profiles_info[code]['XAXIS'][binnum_resolution]

    if plot_kwargs['plot_ind_number_lines']:
        axis.axvline(x=x_resolution, color=color, ls=':', alpha=0.3, zorder=plot_kwargs['number_line_zorder'])
    return x_resolution

    
def plot_property(profiles_info, code, prop, bin_,  type_='regular', classification_='all',
                  element=None, element_num=None, element_denom=None,
                  axis=plt, avg_color='black', ind_color='black',
                  avg_ls='-', ind_ls='-', avg_lw=3, sm=None, 
                  label=None, extra_label=None,
                  fig_kwargs=None, plot_kwargs=None,
                  current_axis_info=None):
    
    # drawstyle = {'default', 'steps', 'steps-pre', 'steps-mid', 'steps-post'}
    
    if plot_kwargs['show_labels']:
        labels = {
#             'std':'std!'+label,
#             'quantiles':'quantiles!'+label,
#             'mean':'mean!'+label,
#             'median':'median!'+label,
            'avg':'avg!'+label,
            'limits':'limits!'+label,
        }
    else:
        labels = {
#             'std':None,
#             'quantiles':None,
#             'mean':None,
#             'median':None,
            'avg':None,
            'limits':None,
        }
        
    original_prop = prop
    if fig_kwargs['yaxis']['is_log']:
        prop = 'log_' + str(prop)
        
    classification = classification_ + '_profiles'
    
    ## Check if classification exists
    try:
        test_ = list(profiles_info[code]['profiles_stats'][bin_][classification]['avg'].keys())[0]
    except:
        if plot_kwargs['verbose']:
            print(code, ':', classification, 'not available')
            print()
        return
    
    ## Check if property exists
    try:
        test_ = profiles_info[code]['profiles_stats'][bin_][classification]['avg'][prop]
    except:
        if plot_kwargs['verbose']:
            print(code, ':', prop, 'not available')
            print()
        return
    
    
    metal_idx_ = None
    metal_idx_num_ = None
    metal_idx_denom_ = None
    
    if type_.lower() == 'regular':
        try:
            test_ = profiles_info[code]['profiles_stats'][bin_][classification]['avg'][prop]
        except:
            if plot_kwargs['verbose']:
                print(code, classification, ':', prop, 'not available')
                print()
            return
        
        vals_ = {
#             'mean':profiles_info[code]['profiles_stats'][bin_][classification]['mean'][prop],
#             'std':profiles_info[code]['profiles_stats'][bin_][classification]['std'][prop],
#             'median':profiles_info[code]['profiles_stats'][bin_][classification]['median'][prop],
#             'quantiles_16':profiles_info[code]['profiles_stats'][bin_][classification]['quantiles'][prop]['16'],
#             'quantiles_84':profiles_info[code]['profiles_stats'][bin_][classification]['quantiles'][prop]['84'],
            'avg':profiles_info[code]['profiles_stats'][bin_][classification]['avg'][prop],
            'limit_lo':profiles_info[code]['profiles_stats'][bin_][classification]['limits'][prop]['lo'],
            'limit_hi':profiles_info[code]['profiles_stats'][bin_][classification]['limits'][prop]['hi'],
        }
    
        
    elif type_.lower() == 'metallicity':
        try:
            metal_idx_ = profiles_info[code]['metals_idx'][element]#-1
        except:
            if plot_kwargs['verbose']:
                print("%s not available for %s" % (element, code))
                print()
            return
        
        vals_ = {
#             'mean':profiles_info[code]['profiles_stats'][bin_][classification]['mean'][prop][:,metal_idx_],
#             'std':profiles_info[code]['profiles_stats'][bin_][classification]['std'][prop][:,metal_idx_],
#             'median':profiles_info[code]['profiles_stats'][bin_][classification]['median'][prop][:,metal_idx_],
#             'quantiles_16':profiles_info[code]['profiles_stats'][bin_][classification]['quantiles'][prop]['16'][:,metal_idx_],
#             'quantiles_84':profiles_info[code]['profiles_stats'][bin_][classification]['quantiles'][prop]['84'][:,metal_idx_],
            'avg':profiles_info[code]['profiles_stats'][bin_][classification]['avg'][prop][:,metal_idx_],
            'limit_lo':profiles_info[code]['profiles_stats'][bin_][classification]['limits'][prop]['lo'][:,metal_idx_],
            'limit_hi':profiles_info[code]['profiles_stats'][bin_][classification]['limits'][prop]['hi'][:,metal_idx_],
        }
        
    elif type_.lower() == 'metallicity ratio':
        try:
            metal_idx_num_ = profiles_info[code]['metals_idx'][element_num]
            metal_idx_denom_ = profiles_info[code]['metals_idx'][element_denom]
        except:
            if plot_kwargs['verbose']:
                print("%s/%s not available for %s" % (element_num, element_denom, code))
                print()
            return
        
        vals_ = {
#             'mean':profiles_info[code]['profiles_stats'][bin_][classification]['mean'][prop][:,metal_idx_num_,metal_idx_denom_],
#             'std':profiles_info[code]['profiles_stats'][bin_][classification]['std'][prop][:,metal_idx_num_,metal_idx_denom_],
#             'median':profiles_info[code]['profiles_stats'][bin_][classification]['median'][prop][:,metal_idx_num_,metal_idx_denom_],
#             'quantiles_16':profiles_info[code]['profiles_stats'][bin_][classification]['quantiles'][prop]['16'][:,metal_idx_num_,metal_idx_denom_],
#             'quantiles_84':profiles_info[code]['profiles_stats'][bin_][classification]['quantiles'][prop]['84'][:,metal_idx_num_,metal_idx_denom_],
            'avg':profiles_info[code]['profiles_stats'][bin_][classification]['avg'][prop][:,metal_idx_num_,metal_idx_denom_],
            'limit_lo':profiles_info[code]['profiles_stats'][bin_][classification]['limits'][prop]['lo'][:,metal_idx_num_,metal_idx_denom_],
            'limit_hi':profiles_info[code]['profiles_stats'][bin_][classification]['limits'][prop]['hi'][:,metal_idx_num_,metal_idx_denom_],
            
        }
        
    else:
        if plot_kwargs['verbose']:
            print('type_ must be regular, metallicity, or metallicity ratio')
        return

    
    
    ## Get correct units
    try:
        og_prop_yunits = vals_[list(vals_.keys())[0]].units
    except:
        og_prop_yunits = profiles_info[code]['profiles'][bin_][classification][profiles_info[code]['profiles'][bin_][classification_+'_id'][0]]['y'][prop].units
        for value in vals_.keys():
            vals_[value] = pnb.array.SimArray(vals_[value], units=prop_yunits)
    
    if plot_kwargs['verbose']:
        print('Original property yunits:', og_prop_yunits)
        print()
        
    if fig_kwargs['yaxis']['use_units']:
        for value in vals_.keys():
            vals_[value] = vals_[value].in_units(fig_kwargs['yaxis']['units'][original_prop])
        prop_yunits = vals_[list(vals_.keys())[0]].units
        if plot_kwargs['verbose']:
            print('New property yunits:', prop_yunits)
    else:
        prop_yunits = og_prop_yunits
    
    
    ## individual profiles
    if plot_kwargs['ind']['plot']:
        plot_individual_profiles(profiles_info, code, prop, bin_, type_=type_, classification_=classification_,
                                 axis=axis, ind_color=ind_color, ind_ls=ind_ls, sm=sm,
                                 metal_idx=metal_idx_, metal_idx_num=metal_idx_num_, metal_idx_denom=metal_idx_denom_,
                                 og_prop_yunits=og_prop_yunits, prop_yunits=prop_yunits,
                                 label=label,
                                 plot_kwargs=plot_kwargs['ind'], fig_kwargs=fig_kwargs, 
                                 current_axis_info=current_axis_info)                 
                
            
    
    ## Profile statistics
    if plot_kwargs['limits']['plot']:
        axis.fill_between(profiles_info[code]['XAXIS'],
                          vals_['limit_lo'], vals_['limit_hi'],
                          color=avg_color, alpha=plot_kwargs['limits']['alpha'], 
                          label=labels['limits'], zorder=plot_kwargs['limits']['zorder'])
        
    if plot_kwargs['avg']['plot']:
        axis.plot(profiles_info[code]['XAXIS'], vals_['avg'], 
                  color=avg_color, ls=avg_ls, lw=avg_lw, label=labels['avg'], #plot_kwargs['avg']['lw'], 
                  drawstyle=plot_kwargs['avg']['drawstyle'], alpha=plot_kwargs['avg']['alpha'], 
                  zorder=plot_kwargs['avg']['zorder'])
    
#     if plot_kwargs['std']['plot']:
#         axis.fill_between(profiles_info[code]['XAXIS'],
#                           vals_['mean'] - vals_['std'],
#                           vals_['mean'] + vals_['std'],
#                           color=avg_color, alpha=plot_kwargs['std']['alpha'], 
#                           label=labels['std'], zorder=plot_kwargs['std']['zorder'])
    
#     if plot_kwargs['quantiles']['plot']:
#         axis.fill_between(profiles_info[code]['XAXIS'],
#                           vals_['quantiles_16'], vals_['quantiles_84'],
#                           color=avg_color, alpha=plot_kwargs['quantiles']['alpha'], 
#                           label=labels['quantiles'], zorder=plot_kwargs['quantiles']['zorder'])

#     if plot_kwargs['mean']['plot']:
#         axis.plot(profiles_info[code]['XAXIS'], vals_['mean'], 
#                   color=avg_color, ls=avg_ls, lw=plot_kwargs['mean']['lw'], label=labels['mean'], 
#                   drawstyle=plot_kwargs['mean']['drawstyle'], alpha=plot_kwargs['mean']['alpha'], 
#                   zorder=plot_kwargs['mean']['zorder'])
        
#     if plot_kwargs['median']['plot']:
#         axis.plot(profiles_info[code]['XAXIS'], vals_['median'], 
#                   color=avg_color, ls=avg_ls, lw=plot_kwargs['median']['lw'], label=labels['median'], 
#                   drawstyle=plot_kwargs['median']['drawstyle'], alpha=plot_kwargs['median']['alpha'], 
#                   zorder=plot_kwargs['median']['zorder'])

    
    if plot_kwargs['show_extra_labels']:
        for extra_label_ in extra_label:
            axis.plot([], [], color=avg_color, ls=avg_ls, lw=plot_kwargs['legend_avg']['lw'], 
                      drawstyle=plot_kwargs['legend_avg']['drawstyle'], alpha=plot_kwargs['legend_avg']['alpha'], 
                      label=extra_label_)




            
            

        
def create_bin_label(bin_name, bin_edges, bin_, log_bin=False, set_num_decimals=False, num_decimals=1, exact=True):
    if exact:
        inequality1 = '<'
        inequality2 = r'$\leq$'
    else:
        inequality1 = r'$\lesssim$'
        inequality2 = r'$\lesssim$'
    
    if log_bin:
        if set_num_decimals:
            label = r'%s%s$\log\,$%s%s%s' % (round(np.log10(bin_edges[bin_]), num_decimals), 
                                             inequality1,
                                             bin_name, 
                                             inequality2,
                                             round(np.log10(bin_edges[bin_+1]), num_decimals))
        else:
            label = r'%s%s$\log\,$%s%s%s' % (np.log10(bin_edges[bin_]), 
                                             inequality1,
                                             bin_name, 
                                             inequality2,
                                             np.log10(bin_edges[bin_+1]))
    else:
        if set_num_decimals:
            label = r'%s%s%s%s%s' % (round(bin_edges[bin_], num_decimals), 
                                     inequality1,
                                     bin_name, 
                                     inequality2,
                                     round(bin_edges[bin_+1], num_decimals))
        else:
            label = r'%s%s%s%s%s' % (bin_edges[bin_],
                                     inequality1,
                                     bin_name,inequality2,
                                     bin_edges[bin_+1])
        
    return label

        
        
        
        
        
def set_grid_val(grids, val, which='overlay'):
    if grids[which]['type'] == 'bins':
        bin_ = val
    elif grids[which]['type'] == 'sims':
        code = val
        info = profiles_info[code]
    elif grids[which]['type'] == 'props':
        prop = val
    elif grids[which]['type'] == 'classifications':
        classification = val
        
    

    


def recursive_loop(total_iter, current_iter=1, count=0,
                   num_per_loop=None, current_result='',
                   func=None, trigger_func=None, 
                   args=[], kwargs={}, 
                   trigger_args=[], trigger_kwargs={}):
    
    if current_iter < total_iter:
        for ii in range(num_per_loop[current_iter-1]):
            result = func(ii, current_iter-1, *args, **kwargs)
            current_result += result + '; '
            
            current_result, count = recursive_loop(total_iter, current_iter=current_iter+1, count=count,
                           num_per_loop=num_per_loop, current_result=current_result,
                           func=func, trigger_func=trigger_func,
                           args=args, kwargs=kwargs,
                           trigger_args=trigger_args, trigger_kwargs=trigger_kwargs)
        
            ## removes last iteration from result
            result_parts = current_result.split('; ')
            current_result = ''
            for part in result_parts[:-2]:
                current_result += part + '; '
            
    elif current_iter == total_iter:
        for jj in range(num_per_loop[current_iter-1]):
            result = func(jj, current_iter-1, *args, **kwargs)
            current_result += result
            
            trigger_func(count, *trigger_args, label=current_result, **trigger_kwargs)

            ## removes last iteration from result
            result_parts = current_result.split('; ')
            current_result = ''
            for part in result_parts[:-1]:
                current_result += part + '; '
                
            count += 1

    return current_result, count
    
    
    
    
    
    
def update_info(iteration, level, current_info, 
                 grids=None, all_grid_types=None,
                 profiles_info=None, props=None, type_='regular',
                 elements=None, elements_num=None, elements_denom=None,
                 plot_kwargs=None, current_axis_info=None):
    
    if not isinstance(level, str):
        level = 'overlay' + str(level)
    
    current_grid = grids[level]
    
    
    ## Get true iteration (ie. if using next_col_after_x_rows or next_row_after_x_cols)
    ## e.g. if fig is 9 rows by 1 column, but next_col_after_x_rows=3 so that it's actually a 3x3 figure,
    ## and we are on the second column with row iteration=2, then the "true" row iteration is actually 5
    if (level == 'row' or level == 'col') and (plot_kwargs['figure']['reshape_width'] or plot_kwargs['figure']['reshape_height']):
        if level == 'row':
            other_level = 'col'
        else:
            other_level = 'row'

        other_iteration = current_axis_info[other_level]
#             print('other_iteration:', other_iteration)
        level_length = grids['n'+level+'s']
#             print('level_length:', level_length)
        true_iteration = iteration + other_iteration * level_length
    
    else:
        true_iteration = copy.deepcopy(iteration)


    
    
    grid_val = current_grid['grid'][true_iteration]
#     if plot_kwargs['display']['verbose']:
#         print(level, true_iteration, ':', current_grid['type'], grid_val)
        
        
    if current_grid['type'] == 'bins':
        bin_ = grid_val
        current_info['bin_'] = bin_

    elif current_grid['type'] == 'sims':
        code = grid_val
        info = profiles_info[code]
        current_info['code'] = code
        current_info['info'] = info
        
    elif current_grid['type'] == 'props':
        prop = grid_val
        current_info['prop'] = prop
        prop_specific_name = plot_kwargs['prop_specific_names'][true_iteration]
        if type_.lower() == 'metallicity':
            if 'elements' not in all_grid_types:
                element = elements[true_iteration]
                current_info['element'] = element
        elif type_.lower() == 'metallicity ratio':
            if 'elements_num' in all_grid_types:
                element_denom = elements_denom[true_iteration]
                current_info['element_denom'] = element_denom
            elif 'elements_denom' in all_grid_types:
                element_num = elements_num[true_iteration]
                current_info['element_num'] = element_num
            else:
                element_num = elements_num[true_iteration]
                element_denom = elements_denom[true_iteration]
                current_info['element_num'] = element_num
                current_info['element_denom'] = element_denom
                
    elif current_grid['type'] == 'classifications':
        classification = grid_val
        current_info['classification'] = classification
        
    elif current_grid['type'] == 'elements':
        element = grid_val
        current_info['element'] = element
        
    elif current_grid['type'] == 'elements_num':
        element_num = grid_val
        current_info['element_num'] = element_num
        
    elif current_grid['type'] == 'elements_denom':
        element_denom = grid_val
        current_info['element_denom'] = element_denom
        
    elif current_grid['type'] == 'elements_num_denom':
        element_num = grid_val[0]
        element_denom = grid_val[1]
        current_info['element_num'] = element_num
        current_info['element_denom'] = element_denom
        
        
        
#     if plot_kwargs['display']['verbose']:
#         print('current_info:')
#         print(current_info.keys())
#         print('code:', current_info['code'])
#         print('bin_:', current_info['bin_'])
#         print('prop:', current_info['prop'])
#         print('classification:', current_info['classification'])
#         print('element:', current_info['element'])
#         print('element_num:', current_info['element_num'])
#         print('element_denom:', current_info['element_denom'])
#         print()
    
    
    
    
    
def update_info_and_create_label(iteration, level, current_info, 
                                 grids=None, all_grid_types=None, bin_info=None,
                                 profiles_info=None, props=None, type_='regular',
                                 elements=None, elements_num=None, elements_denom=None,
                                 plot_kwargs=None, sim_plot_kwargs=None,
                                 current_axis_info=None):
    
    if not isinstance(level, str):
        level = 'overlay' + str(level)
    
    current_grid = grids[level]
    
    
    ## Get true iteration (ie. if using next_col_after_x_rows or next_row_after_x_cols)
    ## e.g. if fig is 9 rows by 1 column, but next_col_after_x_rows=3 so that it's actually a 3x3 figure,
    ## and we are on the second column with row iteration=2, then the "true" row iteration is actually 5
    if level == 'row' or level == 'col':
        if plot_kwargs['figure']['reshape_width']:
    #         if level == 'row' or level == 'col':
            if level == 'row':
                print('row')
                other_level = 'col'
                other_iteration = current_axis_info[other_level]
                print('other_iteration:', other_iteration)
                level_length = grids['n'+level+'s']
                print('level_length:', level_length)
                true_iteration = iteration + other_iteration * level_length
            else:
                print('col')
    #             other_level = 'row'
#                 if recalculate:
                update_info(current_axis_info['row'], 'row', current_info, 
                             grids=grids, all_grid_types=all_grid_types,
                             profiles_info=profiles_info, props=props, type_=type_,
                             elements=elements, elements_num=elements_num, elements_denom=elements_denom,
                             plot_kwargs=plot_kwargs, current_axis_info=current_axis_info)  ## Gonna need to do this either way (row or col), so can put out of if-else clause?
                true_iteration = current_grid['length'] - 1
        
        elif plot_kwargs['figure']['reshape_height']:
            if level == 'col':
                print('col')
                other_level = 'row'
                other_iteration = current_axis_info[other_level]
                print('other_iteration:', other_iteration)
                level_length = grids['n'+level+'s']
                print('level_length:', level_length)
                true_iteration = iteration + other_iteration * level_length
            else:
                print('row')
    #             other_level = 'row'
#                 if recalculate:
                update_info(current_axis_info['col'], 'col', current_info, 
                             grids=grids, all_grid_types=all_grid_types,
                             profiles_info=profiles_info, props=props, type_=type_,
                             elements=elements, elements_num=elements_num, elements_denom=elements_denom,
                             plot_kwargs=plot_kwargs, current_axis_info=current_axis_info)  ## Gonna need to do this either way (row or col), so can put out of if-else clause?
                true_iteration = current_grid['length'] - 1
        
        
#         elif plot_kwargs['figure']['reshape_height']:
#             if level == 'col':
#                 other_level = 'row'
#                 other_iteration = current_axis_info[other_level]
#                 level_length = grids['n'+level+'s']
#                 true_iteration = iteration + other_iteration * level_length
#             else:
#     #             other_level = 'col'
#                 true_iteration = current_grid['length'] - 1
        
        else:
            true_iteration = copy.deepcopy(iteration)
            
    else:
        true_iteration = copy.deepcopy(iteration)

    if plot_kwargs['display']['verbose']:
        print(level, true_iteration, ':', current_grid['type'])


    
    
    grid_val = current_grid['grid'][true_iteration]
    if plot_kwargs['display']['verbose']:
        print(level, true_iteration, ':', current_grid['type'], grid_val)
        
    prefix = ''
        
    label = None
    if current_grid['type'] == 'bins':
        bin_ = grid_val
        current_info['bin_'] = bin_
        label = prefix + create_bin_label(bin_info['bin_name'], bin_info['bin_edges'], bin_, 
                                                    log_bin=plot_kwargs['binning']['log_bin'],
                                                    set_num_decimals=plot_kwargs['binning']['set_num_decimals'],
                                                    num_decimals=plot_kwargs['binning']['num_decimals'])

        if plot_kwargs['binning']['show_other_bin']:
            other_label = create_bin_label(plot_kwargs['binning']['other_bin_name'], 
                                                   bin_info['other_bin_edges'], bin_, 
                                                   log_bin=plot_kwargs['binning']['log_other_bin'],
                                                   set_num_decimals=plot_kwargs['binning']['set_num_decimals_other_bin'],
                                                   num_decimals=plot_kwargs['binning']['num_decimals_other_bin'], 
                                                   exact=False)
            label += '\n' + '(' + other_label + ')'

    elif current_grid['type'] == 'sims':
        code = grid_val
        info = profiles_info[code]
        current_info['code'] = code
        current_info['info'] = info
        label = prefix + code
        
    elif current_grid['type'] == 'props':
        prop = grid_val
        current_info['prop'] = prop
        prop_specific_name = plot_kwargs['prop_specific_names'][true_iteration]
        if type_.lower() == 'regular':
            label = prefix + prop_specific_name
#             if 'elements' not in all_grid_types:
#                 element = elements[iteration]
#                 current_info['element'] = element
#                 if plot_kwargs['grid']['show_full_prop_name']:
#                     label = prefix + prop_specific_name + ' [%s]' % element
#                 else:
#                     label = prefix + element
#             else:
#                 label = prefix + prop_specific_name
        elif type_.lower() == 'metallicity':
            if 'elements' not in all_grid_types:
                element = elements[true_iteration]
                current_info['element'] = element
                if plot_kwargs['grid']['show_full_prop_name']:
                    label = prefix + prop_specific_name + ' [%s]' % element
                else:
                    label = prefix + element
            else:
                label = prefix + prop_specific_name
        elif type_.lower() == 'metallicity ratio':
            if ('elements_num' in all_grid_types and 'elements_denom' in all_grid_types) or ('elements_num_denom' in all_grid_types):
                label = prefix + prop_specific_name
            elif 'elements_num' in all_grid_types:
                element_denom = elements_denom[true_iteration]
                current_info['element_denom'] = element_denom
                if plot_kwargs['grid']['show_full_prop_name']:
                    label = prefix + prop_specific_name + ' [%s]' % element_denom
                elif plot_kwargs['grid']['show_name'].lower() == 'element':
                    label = prefix + element_denom
                elif plot_kwargs['grid']['show_name'].lower() == 'prop':
                    label = prefix + prop_specific_name
            elif 'elements_denom' in all_grid_types:
                element_num = elements_num[true_iteration]
                current_info['element_num'] = element_num
                if plot_kwargs['grid']['show_full_prop_name']:
                    label = prefix + prop_specific_name + ' [%s]' % element_num
                elif plot_kwargs['grid']['show_name'].lower() == 'element':
                    label = prefix + element_num
                elif plot_kwargs['grid']['show_name'].lower() == 'prop':
                    label = prefix + prop_specific_name
            else:
                element_num = elements_num[true_iteration]
                element_denom = elements_denom[true_iteration]
                current_info['element_num'] = element_num
                current_info['element_denom'] = element_denom
                if plot_kwargs['grid']['show_full_prop_name']:
                    label = prefix + prop_specific_name + ' [%s/%s]' % (element_num, element_denom)
                elif plot_kwargs['grid']['show_name'].lower() == 'element':
                    label = prefix + element_num + '/' + element_denom
                elif plot_kwargs['grid']['show_name'].lower() == 'prop':
                    label = prefix + prop_specific_name
                
            
#             if 'elements_num' not in all_grid_types and 'elements_denom' not in all_grid_types:
#                 element_num = elements_num[iteration]
#                 element_denom = elements_denom[iteration]
#                 current_info['element_num'] = element_num
#                 current_info['element_denom'] = element_denom
#                 if plot_kwargs['grid']['show_full_prop_name']:
#                     label = prefix + prop_specific_name + ' [%s/%s]' % (element_num, element_denom)
#                 else:
#                     label = prefix + element_num + '/' + element_denom
#             elif 'elements_num' not in all_grid_types:
#                 element_num = elements_num[iteration]
#                 current_info['element_num'] = element_num
#                 if plot_kwargs['grid']['show_full_prop_name']:
#                     label = prefix + prop_specific_name + ' [%s]' % element_num
#                 else:
#                     label = prefix + element_num
#             elif 'elements_denom' not in all_grid_types:
#                 element_denom = elements_denom[iteration]
#                 current_info['element_denom'] = element_denom
#                 if plot_kwargs['grid']['show_full_prop_name']:
#                     label = prefix + prop_specific_name + ' [%s]' % element_denom
#                 else:
#                     label = prefix + element_denom
#             else:
#                 label = prefix + prop_specific_name
                
    elif current_grid['type'] == 'classifications':
        classification = grid_val
        current_info['classification'] = classification
        label = prefix + sim_plot_kwargs['classifications']['labels'][true_iteration]
        
    elif current_grid['type'] == 'elements':
        element = grid_val
        current_info['element'] = element
        label = prefix + element
    elif current_grid['type'] == 'elements_num':
        element_num = grid_val
        current_info['element_num'] = element_num
        label = prefix + element_num
    elif current_grid['type'] == 'elements_denom':
        element_denom = grid_val
        current_info['element_denom'] = element_denom
        label = prefix + element_denom
    elif current_grid['type'] == 'elements_num_denom':
        element_num = grid_val[0]
        element_denom = grid_val[1]
        current_info['element_num'] = element_num
        current_info['element_denom'] = element_denom
        label = prefix + element_num + '/' + element_denom
        
        
        
    if plot_kwargs['display']['verbose']:
        print('current_info:')
        print(current_info.keys())
        print('code:', current_info['code'])
        print('bin_:', current_info['bin_'])
        print('prop:', current_info['prop'])
        print('classification:', current_info['classification'])
        print('element:', current_info['element'])
        print('element_num:', current_info['element_num'])
        print('element_denom:', current_info['element_denom'])
        print()
        
        
    return label


        
        
        
def plot_overlay(overlay_tot, current_info, current_axis_info=None, 
                 label=None, ax=plt, sm=None,
                 profiles_info=None, type_='regular',
                 plot_kwargs=None, sim_plot_kwargs=None):
    
#     if not sim_plot_kwargs['really_plot?'][overlay_tot]:
    if not sim_plot_kwargs['really_plot?'][current_axis_info['row']][current_axis_info['col']][overlay_tot]:
        if sim_plot_kwargs['verbose']:
            print('Not plotting overlay', overlay_tot)
        return
    
    if current_info['bin_'] not in current_info['info']['profiles']['bins']:
        if plot_kwargs['display']['verbose']:
            print('bin', current_info['bin_'], 'not available in', current_info['code'])
        return
    
#     label = '#sim##' + label
    label = profiles_info[current_info['code']]['LABEL'] + label
    
    color = sim_plot_kwargs['color_list'][overlay_tot]
    avg_line_style = sim_plot_kwargs['legend_avg']['line_styles'][overlay_tot]
    ind_line_style = sim_plot_kwargs['ind']['line_styles'][overlay_tot]
    avg_lw = sim_plot_kwargs['avg']['lw'][overlay_tot]

    avg_color = color
    ind_color = 'black'
    if sim_plot_kwargs['classifications']['use_colors']:
        classification_color = sim_plot_kwargs['classifications']['colors'][overlay_tot]
        avg_color = classification_color
        ind_color = classification_color
    if sim_plot_kwargs['ind']['color_same_as_avg']:
        ind_color = avg_color

    extra_label = [label]


    ## Just send current_info, not all the separate things
    plot_property(profiles_info, current_info['code'], current_info['prop'], current_info['bin_'], 
                  type_=type_, classification_=current_info['classification'],
                  element=current_info['element'], element_num=current_info['element_num'], element_denom=current_info['element_denom'], 
                  axis=ax, avg_color=avg_color, ind_color=ind_color,
                  avg_ls=avg_line_style, ind_ls=ind_line_style, avg_lw=avg_lw, sm=sm, 
                  label=label, extra_label=extra_label,
                  fig_kwargs=plot_kwargs, plot_kwargs=sim_plot_kwargs, 
                  current_axis_info=current_axis_info)
    
    return



# def plot_overlay_residual(overlay_tot, current_info, current_axis_info=None, 
#                        label=None, ax=plt, sm=None,
#                        profiles_info=None, type_='regular',
#                        plot_kwargs=None, sim_plot_kwargs=None):
    
#     if not sim_plot_kwargs['really_plot?'][overlay_tot]:
#         if sim_plot_kwargs['verbose']:
#             print('Not plotting overlay', overlay_tot)
#         return
    
#     if current_info['bin_'] not in current_info['info']['profiles']['bins']:
#         if plot_kwargs['display']['verbose']:
#             print('bin', current_info['bin_'], 'not available in', current_info['code'])
#         return
    
#     label = '#sim##' + label
    
#     color = sim_plot_kwargs['residuals']['color_list'][overlay_tot]
#     avg_line_style = sim_plot_kwargs['residuals']['ls'][overlay_tot]
#     ind_line_style = sim_plot_kwargs['ind']['line_styles'][overlay_tot]
#     avg_lw = sim_plot_kwargs['residuals']['lw'][overlay_tot]

#     avg_color = color
#     ind_color = 'black'
#     if sim_plot_kwargs['classifications']['use_colors']:
#         classification_color = sim_plot_kwargs['classifications']['colors'][overlay_tot]
#         avg_color = classification_color
#         ind_color = classification_color
#     if sim_plot_kwargs['ind']['color_same_as_avg']:
#         ind_color = avg_color

#     extra_label = [label]


#     ## Just send current_info, not all the separate things
#     plot_property_residual(profiles_info, current_info['code'], current_info['prop'], current_info['bin_'], 
#                   type_=type_, classification_=current_info['classification'],
#                   element=current_info['element'], element_num=current_info['element_num'], element_denom=current_info['element_denom'], 
#                   axis=ax, avg_color=avg_color, ind_color=ind_color,
#                   avg_ls=avg_line_style, ind_ls=ind_line_style, avg_lw=avg_lw, sm=sm, 
#                   label=label, extra_label=extra_label,
#                   fig_kwargs=plot_kwargs, plot_kwargs=sim_plot_kwargs, 
#                   current_axis_info=current_axis_info)
    
#     return
    
    
    
        
        
    
# def calc_ratio(a, b):
#     return a/b

# def calc_frac_residual(a, b):
#     return (a-b)/b

# residuals_dict = {
#     'ratio':calc_ratio,
#     'frac_residual':calc_frac_residual,
# }
        
        
        
        
        
        
        
## Plot in a grid
def plot_grid(profiles_info, props, type_='regular', #units=None,
              elements=None, elements_num=None, elements_denom=None, 
              plot_kwargs=None, sim_plot_kwargs=None, obs_plot_kwargs=None):    
    
#     ## Use your discretion to make sure each simulation has its bin numbers aligning with the same bin
#     ## edges as all the other simulations!!

    bin_edges = profiles_info[list(profiles_info.keys())[0]]['BIN_EDGES']
    bin_name = profiles_info[list(profiles_info.keys())[0]]['BIN_NAME']

    bin_info = {
        'bin_edges':bin_edges,
        'bin_name':bin_name,
    }

    if plot_kwargs['binning']['show_other_bin']:
        if plot_kwargs['display']['verbose']:
            print('Showing alternate bin prop:', plot_kwargs['binning']['other_prop'])
            
        try:
            other_bin_edges = plot_kwargs['binning']['other_bin_func'](bin_edges)
            bin_info['other_bin_edges'] = other_bin_edges
        except BaseException:
            logging.exception('Alternate bin prop could not be calculated correctly')
    
    
    grids = {
        'row':{
            'type':plot_kwargs['grid']['row'],
            'grid':None,
            'length':0,
        },
        'col':{
            'type':plot_kwargs['grid']['col'],
            'grid':None,
            'length':0,
        },
    }
    
    for ii in range(len(plot_kwargs['grid']['overlay'])):
        grids['overlay' + str(ii)] = {
            'type':plot_kwargs['grid']['overlay'][ii],
            'grid':None,
            'length':0,
        }
    
    
    
    all_grid_types = []
    for key, vals in grids.items():
        all_grid_types.append(vals['type'])
        if plot_kwargs['display']['verbose']:
            print(key)
            print(vals['type'])
        if vals['type'] == 'bins':
            vals['grid'] = sim_plot_kwargs['bins_to_plot']
        elif vals['type'] == 'sims':
            vals['grid'] = list(profiles_info.keys())
        elif vals['type'] == 'props':
            vals['grid'] = props
        elif vals['type'] == 'classifications':
            vals['grid'] = sim_plot_kwargs['classifications']['which']
        elif vals['type'] == 'elements':
            vals['grid'] = elements
        elif vals['type'] == 'elements_num':
            vals['grid'] = elements_num
        elif vals['type'] == 'elements_denom':
            vals['grid'] = elements_denom
        elif vals['type'] == 'elements_num_denom':
#             vals['grid'] = [elements_num, elements_denom]
            vals['grid'] = [[element_num, element_denom] for element_num, element_denom in zip(elements_num, elements_denom)]
        else:
            raise Exception("%s must be 'bins', 'sims', 'props', or 'classifications'" % (key+'_'+'grid_type'))
        if plot_kwargs['display']['verbose']:
            print(vals['grid'])
            print()
        vals['length'] = len(vals['grid'])
    
    
    
    ## determine nrows and ncols for figure
    grids['nrows'] = grids['row']['length']
    grids['ncols'] = grids['col']['length']
    
    ## This is only for if the original figure size is 1xN or Nx1, to reshape it into IxJ grid
    ## SO: make sure that "second" dimension of figure has next...after_x... bigger than col length or row length
    if plot_kwargs['figure']['reshape_width']:
        grids['nrows'] = min(plot_kwargs['figure']['next_col_after_x_rows'], grids['row']['length'])
        grids['ncols'] = int(np.ceil(grids['row']['length']/plot_kwargs['figure']['next_col_after_x_rows']))
    
    if plot_kwargs['figure']['reshape_height']:
        grids['ncols'] = min(plot_kwargs['figure']['next_row_after_x_cols'], grids['col']['length'])
        grids['nrows'] = int(np.ceil(grids['col']['length']/plot_kwargs['figure']['next_row_after_x_cols']))
    
    
    
    if plot_kwargs['display']['verbose']:
        print()
        print(grids)
        print()
        
        
        
        
        

    code = list(profiles_info.keys())[sim_plot_kwargs['unused']['sim']]
    info = profiles_info[code]
    bin_ = sim_plot_kwargs['bins_to_plot'][sim_plot_kwargs['unused']['bin']]
    prop = props[sim_plot_kwargs['unused']['prop']]
    prop_specific_name = plot_kwargs['prop_specific_names'][sim_plot_kwargs['unused']['prop']]
    classification = sim_plot_kwargs['classifications']['which'][sim_plot_kwargs['unused']['classification']]
    
    if type_.lower() == 'regular':
        element = None
        element_num = None
        element_denom = None
    if type_.lower() == 'metallicity':
        if 'elements' not in all_grid_types:
            element = elements[sim_plot_kwargs['unused']['prop']]
        else:
            element = elements[sim_plot_kwargs['unused']['element']]
        element_num = None
        element_denom = None
    if type_.lower() == 'metallicity ratio':
        element = None
        if ('elements_num' in all_grid_types and 'elements_denom' in all_grid_types) or ('elements_num_denom' in all_grid_types):
            element_num = elements_num[sim_plot_kwargs['unused']['element_num']]
            element_denom = elements_denom[sim_plot_kwargs['unused']['element_denom']]
#             element_num = elements_num[sim_plot_kwargs['unused']['prop']]
#             element_denom = elements_denom[sim_plot_kwargs['unused']['prop']]
        elif 'elements_num' in all_grid_types:
            element_num = elements_num[sim_plot_kwargs['unused']['element_num']]
            element_denom = elements_denom[sim_plot_kwargs['unused']['prop']]
        elif 'elements_denom' in all_grid_types:
            element_num = elements_num[sim_plot_kwargs['unused']['prop']]
            element_denom = elements_denom[sim_plot_kwargs['unused']['element_denom']]
        else:
            element_num = elements_num[sim_plot_kwargs['unused']['prop']]
            element_denom = elements_denom[sim_plot_kwargs['unused']['prop']]
            
            
            
    current_info = {
        'code':code,
        'info':info,
        'bin_':bin_,
        'prop':prop,
        'classification':classification,
        'element':element,
        'element_num':element_num,
        'element_denom':element_denom,
    }
    
    current_axis_info = {
        'row':0,
        'col':0,
        'overlay_tot':0,
    }
    
    for ii in range(len(plot_kwargs['grid']['overlay'])):
        current_axis_info['overlay'+str(ii)] = 0
    
    
    
    if plot_kwargs['cbar']['plot']:
        cbarlabel = bin_name
        if plot_kwargs['cbar']['log_cbar']:
            cbarlabel = r'$\log \,$%s' % cbarlabel


#     ## Add dummy subplots to right or bottom of figure
#     ## Allows to loop over both rows and columns regardless of the dimensionality of the figure
#     edge_axes_to_remove = []
#     if plot_kwargs['display']['extra_empty_row_subplot'] or plot_kwargs['display']['extra_empty_col_subplot']:
#         if plot_kwargs['display']['extra_empty_row_subplot'] and not plot_kwargs['display']['extra_empty_col_subplot']:
#             gridspec_kw={
#                 'wspace':plot_kwargs['figure']['col_spacing'], 
#                 'hspace':plot_kwargs['figure']['row_spacing'], 
#                 'height_ratios':list([1]*grids['row']['length'])+[1e-30],
#             }
#             fig, ax = plt.subplots(nrows=grids['row']['length']+1, ncols=grids['col']['length'], 
#                                    sharex=plot_kwargs['figure']['sharex'], sharey=plot_kwargs['figure']['sharey'],
#                                    figsize=(plot_kwargs['figure']['width'], plot_kwargs['figure']['height']), 
#                                    constrained_layout=plot_kwargs['figure']['constrained_layout'], gridspec_kw=gridspec_kw,
#                                    num=1, clear=True, facecolor=plot_kwargs['figure']['facecolor'])
            
#             edge_axes_to_remove.append('row')
# #             for _ax_ in ax[-1,:]:
# #                 _ax_.axis('off')

# #             nrows = ax.shape[0]-1
# #             ncols = ax.shape[1]
#         elif (not plot_kwargs['display']['extra_empty_row_subplot']) and plot_kwargs['display']['extra_empty_col_subplot']:
#             gridspec_kw={
#                 'wspace':plot_kwargs['figure']['col_spacing'], 
#                 'hspace':plot_kwargs['figure']['row_spacing'], 
#                 'width_ratios':list([1]*grids['col']['length'])+[1e-30],
#             }
#             fig, ax = plt.subplots(nrows=grids['row']['length'], ncols=grids['col']['length']+1, 
#                                    sharex=plot_kwargs['figure']['sharex'], sharey=plot_kwargs['figure']['sharey'],
#                                    figsize=(plot_kwargs['figure']['width'], plot_kwargs['figure']['height']), 
#                                    constrained_layout=plot_kwargs['figure']['constrained_layout'], gridspec_kw=gridspec_kw,
#                                    num=1, clear=True, facecolor=plot_kwargs['figure']['facecolor'])
            
#             edge_axes_to_remove.append('col')
# #             for _ax_ in ax[:,-1]:
# #                 _ax_.axis('off')

# #             nrows = ax.shape[0]
# #             ncols = ax.shape[1]-1
#         elif plot_kwargs['display']['extra_empty_row_subplot'] and plot_kwargs['display']['extra_empty_col_subplot']:
#             gridspec_kw={
#                 'wspace':plot_kwargs['figure']['col_spacing'], 
#                 'hspace':plot_kwargs['figure']['row_spacing'], 
#                 'height_ratios':list([1]*grids['row']['length'])+[1e-30],
#                 'width_ratios':list([1]*grids['col']['length'])+[1e-30],
# #                 'width_ratios':list([1]*grids['col']['length'])+[1e-30],
# #                 'height_ratios':list([1]*grids['row']['length'])+[1e-30],
#             }
#             fig, ax = plt.subplots(nrows=grids['row']['length']+1, ncols=grids['col']['length']+1, 
#                                    sharex=plot_kwargs['figure']['sharex'], sharey=plot_kwargs['figure']['sharey'],
#                                    figsize=(plot_kwargs['figure']['width'], plot_kwargs['figure']['height']), 
#                                    constrained_layout=plot_kwargs['figure']['constrained_layout'], gridspec_kw=gridspec_kw,
#                                    num=1, clear=True, facecolor=plot_kwargs['figure']['facecolor'])
            
#             edge_axes_to_remove.append('row')
#             edge_axes_to_remove.append('col')
# #             for _ax_ in ax[-1,:]:
# #                 _ax_.axis('off')
                
# #             for _ax_ in ax[:,-1]:
# #                 _ax_.axis('off')

# #             nrows = ax.shape[0]-1
# #             ncols = ax.shape[1]-1
#     else:
#         gridspec_kw={
#             'wspace':plot_kwargs['figure']['col_spacing'], 
#             'hspace':plot_kwargs['figure']['row_spacing'], 
#         }
#         fig, ax = plt.subplots(nrows=grids['row']['length'], ncols=grids['col']['length'], 
#                                sharex=plot_kwargs['figure']['sharex'], sharey=plot_kwargs['figure']['sharey'],
#                                figsize=(plot_kwargs['figure']['width'], plot_kwargs['figure']['height']), 
#                                constrained_layout=plot_kwargs['figure']['constrained_layout'], gridspec_kw=gridspec_kw,
#                                num=1, clear=True, facecolor=plot_kwargs['figure']['facecolor'])
        
        
        
    gridspec_kw={
        'wspace':plot_kwargs['figure']['col_spacing'], 
        'hspace':plot_kwargs['figure']['row_spacing'], 
    }
    fig, ax = plt.subplots(nrows=grids['nrows'], ncols=grids['ncols'], 
                           sharex=plot_kwargs['figure']['sharex'], sharey=plot_kwargs['figure']['sharey'],
                           figsize=(plot_kwargs['figure']['width'], plot_kwargs['figure']['height']), 
                           constrained_layout=plot_kwargs['figure']['constrained_layout'], gridspec_kw=gridspec_kw,
                           num=1, clear=True, facecolor=plot_kwargs['figure']['facecolor'])
        
        
        
#         nrows = ax.shape[0]
#         ncols = ax.shape[1]

#     print(ax)
    
#     if np.ndim(ax) == 1:
#         ax = [ax]
        
    while np.ndim(ax) < 2:
#         print('HELLO!')
        ax = np.array([ax])
    nrows = ax.shape[0]
    ncols = ax.shape[1]
        
        
#     ## Turn off displaying of dummy axes
#     if 'row' in edge_axes_to_remove:
#         for _ax_ in ax[-1,:]:
#             _ax_.axis('off')
#         nrows = ax.shape[0]-1
#     else:
#         nrows = ax.shape[0]
        
#     if 'col' in edge_axes_to_remove:
#         for _ax_ in ax[:,-1]:
#             _ax_.axis('off')
#         ncols = ax.shape[1]-1
#     else:
#         ncols = ax.shape[1]
        
#     print(ncols)
#     print(nrows)
#     print(ax)
        
    
    noverlay = [grids[grid_level]['length'] for grid_level in grids.keys() if 'overlay' in grid_level]
    
    
    ## For plotting residuals
    if sim_plot_kwargs['residuals']['plot']:
#         for _code1_ in profiles_info.keys():
#             _xaxis1_ = profiles_info[_code1_]['XAXIS']
#             for _code2_ in profiles_info.keys():
#                 _xaxis2_ = profiles_info[_code2_]['XAXIS']
#                 assert np.array_equal(_xaxis1_, _xaxis2_, equal_nan=True) , 'Overlaid simulation profiles have different x-axes. Cannot calculate residuals.'               
        
        gs_xaxis = profiles_info[current_info['code']]['XAXIS']
        gs_outer = fig.add_gridspec(nrows=nrows, ncols=ncols, **gridspec_kw) ## !!!!
        
        gs_axis_list = []
        gs_axis_dict = {
            'main_ax':np.full((nrows,ncols), None, dtype=object),
            'residual_ax':np.full((nrows,ncols), None, dtype=object),
            'inner_ax':np.full((2*nrows,ncols), None, dtype=object),
            'outer_ax':np.full((nrows,ncols), None, dtype=object),
        }
        
        
#         if 'row' in edge_axes_to_remove:
#             for _ax_ in gs_outer[-1,:]:
#                 _ax_.axis('off')

#         if 'col' in edge_axes_to_remove:
#             for _ax_ in gs_outer[:,-1]:
#                 _ax_.axis('off')
                
#         current_axis_info['gs_row'] = None

    
    if plot_kwargs['display']['verbose']:
        print()
        print('nrows =', nrows)
        print('ncols =', ncols)
        print('noverlay =', noverlay)
        print()
        
    if plot_kwargs['title']['show']:
        fig.suptitle(plot_kwargs['title']['text'], y=plot_kwargs['title']['y'], va='top', 
                     wrap=plot_kwargs['title']['wrap'], fontsize=plot_kwargs['title']['fontsize'], fontweight='bold')
        

#     fig_cbar = False



    if plot_kwargs['cbar']['plot']:
        if plot_kwargs['cbar']['placement'].lower() == 'ind':
            sm_arr = np.full((nrows,ncols), None, dtype=object)
        elif plot_kwargs['cbar']['placement'].lower() == 'row':
            sm_arr = np.full(nrows, None, dtype=object)
        elif plot_kwargs['cbar']['placement'].lower() == 'col':
            sm_arr = np.full(ncols, None, dtype=object)
        elif plot_kwargs['cbar']['placement'].lower() == 'fig':
            sm_arr = None

        

    for row in range(nrows):
        current_axis_info['row'] = row
        
        row_label = update_info_and_create_label(row, 'row', current_info, 
                                                 grids=grids, all_grid_types=all_grid_types, bin_info=bin_info,
                                                 profiles_info=profiles_info, props=props, type_=type_,
                                                 elements=elements, elements_num=elements_num, elements_denom=elements_denom,
                                                 plot_kwargs=plot_kwargs, sim_plot_kwargs=sim_plot_kwargs,
                                                 current_axis_info=current_axis_info)
        
        if sim_plot_kwargs['residuals']['plot']:
            gs_axis_list.append([])

        
        for col in range(ncols):
            current_axis_info['col'] = col
            
            col_label = update_info_and_create_label(col, 'col', current_info, 
                                                     grids=grids, all_grid_types=all_grid_types, bin_info=bin_info,
                                                     profiles_info=profiles_info, props=props, type_=type_,
                                                     elements=elements, elements_num=elements_num, elements_denom=elements_denom,
                                                     plot_kwargs=plot_kwargs, sim_plot_kwargs=sim_plot_kwargs,
                                                     current_axis_info=current_axis_info)
            
            if sim_plot_kwargs['residuals']['plot']:
                gs_axis_list[row].append([])
            
            
                
#             if plot_kwargs['figure']['sharex'] is not False and plot_kwargs['figure']['sharex'] is not 'none':
#                 print('xxxxx')
#                 if row == nrows:
#                     ax[row,col].tick_params(labelbottom=True)
#                 else:
#                     ax[row,col].tick_params(labelbottom=False)
#             elif plot_kwargs['figure']['hide_xticklabels'] and row != nrows-1:
#                 print('HII')
#                 ax[row,col].tick_params(labelbottom=False)
            
            
            
            ## If plotting a subplot of each subplot that shows the ratios or fractional differences between sims
            if sim_plot_kwargs['residuals']['plot']:
                gs_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs_outer[row,col], 
                                                            hspace=0, height_ratios=plot_kwargs['yaxis']['residuals']['height_ratios'])
#                 gs_inner = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=ax[row,col].get_subplotspec(), 
#                                                             hspace=0, height_ratios=plot_kwargs['yaxis']['residuals']['height_ratios'])
    
                if plot_kwargs['figure']['sharey'] == False or plot_kwargs['figure']['sharey'] == 'none':# or col == 0:
#                     print('col == 0')
                    main_ax = fig.add_subplot(gs_inner[0])
                    residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax)
                elif plot_kwargs['figure']['sharey'] == 'row':
#                     print('col != 0')
                    if col == 0:
                        main_ax = fig.add_subplot(gs_inner[0])
                        residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax)
                    else:
                        main_ax = fig.add_subplot(gs_inner[0], sharey=gs_axis_list[row][0]['main_ax'])
                        residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax, sharey=gs_axis_list[row][0]['residual_ax'])
    #                     main_ax.set_yticklabels([])
    #                     residual_ax.set_yticklabels([])
                        main_ax.tick_params(labelleft=False)
                        residual_ax.tick_params(labelleft=False)
                elif plot_kwargs['figure']['sharey'] == 'col':
                    if row == 0:
                        main_ax = fig.add_subplot(gs_inner[0])
                        residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax)
                    else:
                        main_ax = fig.add_subplot(gs_inner[0], sharey=gs_axis_list[0][col]['main_ax'])
                        residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax, sharey=gs_axis_list[row][0]['residual_ax'])
                elif plot_kwargs['figure']['sharey'] == True or plot_kwargs['figure']['sharey'] == 'all':
                    if row == 0 and col == 0:
                        main_ax = fig.add_subplot(gs_inner[0])
                        residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax)
                    else:
                        main_ax = fig.add_subplot(gs_inner[0], sharey=gs_axis_list[0][0]['main_ax'])
                        residual_ax = fig.add_subplot(gs_inner[1], sharex=main_ax, sharey=gs_axis_list[0][0]['residual_ax'])
                        
                        if col != 0:
                            main_ax.tick_params(labelleft=False)
                            residual_ax.tick_params(labelleft=False)
                            
                if (plot_kwargs['figure']['sharex'] != False and 
                    plot_kwargs['figure']['sharex'] != 'none' and 
                    plot_kwargs['figure']['sharex'] != 'row'):
#                     last_row = 
#                     if plot_kwargs['display']['extra_empty_row_subplot']:
                        
                    if row == nrows-1:
#                         print()
#                         print('nrows:',nrows)
#                         print('last row:',row)
#                         print()
                        ax[row,col].tick_params(labelbottom=True)
                        residual_ax.tick_params(labelbottom=True)
                    else:
#                         print()
#                         print('row:',row)
#                         print()
                        ax[row,col].tick_params(labelbottom=False)
                        residual_ax.tick_params(labelbottom=False)
                    
                            
#                 if plot_kwargs['figure']['sharex'] is not False and plot_kwargs['figure']['sharex'] is not 'none':
#                     print('XXXXX')
#                     if row == nrows:
#                         residual_ax.tick_params(labelbottom=True)
#                     else:
#                         residual_ax.tick_params(labelbottom=False)
#                 elif plot_kwargs['figure']['hide_xticklabels'] and rows != nrows:
#                     residual_ax.tick_params(labelbottom=False)
                    
                    
#                 ax[row,col].set_visible(False)
                
                gs_axis_list[row][col] = {'main_ax':main_ax, 'residual_ax':residual_ax}
                gs_axis_dict['main_ax'][row,col] = main_ax
                gs_axis_dict['residual_ax'][row,col] = residual_ax
#                 gs_axis_dict['inner_ax'][row,col] = gs_inner
                gs_axis_dict['inner_ax'][2*row,col] = main_ax
                gs_axis_dict['inner_ax'][2*row+1,col] = residual_ax
                gs_axis_dict['outer_ax'][row,col] = gs_outer[row,col]
    
#                 if plot_kwargs['figure']['sharex'] is not False and plot_kwargs['figure']['sharex'] is not 'none':
#                     sharex = ax[row,col]
#                 else:
#                     sharex = None
#                 if plot_kwargs['figure']['sharey'] is not False or plot_kwargs['figure']['sharey'] is not 'none':
#                     sharey = ax[row,col]
#                 else:
#                     sharey = None
#                 main_ax = fig.add_subplot(gs_inner[0],
#                                           sharex=sharex, 
#                                           sharey=sharey)
#                 residual_ax = fig.add_subplot(gs_inner[1],
#                                               sharex=sharex, 
#                                               sharey=sharey)
            else:
                main_ax = ax[row,col]
                residual_ax = None
                
                if plot_kwargs['figure']['sharex'] != False and plot_kwargs['figure']['sharex'] != 'none' and row == nrows-1:
                    main_ax.tick_params(labelbottom=True)
                    
                    
                if (plot_kwargs['figure']['sharex'] != False and 
                    plot_kwargs['figure']['sharex'] != 'none' and 
                    plot_kwargs['figure']['sharex'] != 'row'):                        
                    if row == nrows-1:
                        main_ax.tick_params(labelbottom=True)
                    else:
                        main_ax.tick_params(labelbottom=False)
            
                
                
                
            ## Should prob go in overlay loop? --> no, second looop --> no, overlay loop
            ## determine min and max binning values (eg. M500, Tspec,corr, ...)
            ## Can just make it so I provide the subplots I want the colorbar to be drawn on, 
            ## and whether it's rows, columns or individual subplots that the colorbar should include the properties of
            if plot_kwargs['cbar']['plot']:
                
                if plot_kwargs['cbar']['provide_sm']:
                    sm = plot_kwargs['cbar']['sm']
                
                else:
                    ## ADD IN NEW OPTIONS FOR GRID (ELEMENTS, ELEMENTS_NUM, ELEMENTS_DENOM)

                    ## Whether to use the defined bin limits or the min
                    ## and max values of the binning property of all halos in this bin
                    if plot_kwargs['cbar']['method']=='default':
                        bin_color_limits = [bin_edges[bin_], bin_edges[bin_+1]]
                    elif plot_kwargs['cbar']['method']=='provide':
                        bin_color_limits = plot_kwargs['cbar']['limits']
                    elif plot_kwargs['cbar']['method']=='calculate':
                        bin_vals_ = []

                        ## calculate_on can be 'ind', 'row', 'col', or 'all'
                        for overlay_loop_idx in range(noverlay):
                            overlay_val_ = grids['overlay']['grid'][overlay_loop_idx]
                            set_grid_val(grids, overlay_val_, which='overlay')

                            if plot_kwargs['cbar']['calculate_on'] == 'ind':
                                for halo in info['profiles'][bin_][classification+'_id']:
                                    halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                    bin_val_ = info_['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                    bin_vals_.append(bin_val_)


                            elif plot_kwargs['cbar']['calculate_on'] == 'row' or plot_kwargs['cbar']['calculate_on'] == 'col':
                                ## run over all profile binning values in a row or column
                                if grids[plot_kwargs['cbar']['calculate_on']]['type'] == 'sims':
                                    for code, info in profiles_info.items():
                                        for halo in info['profiles'][bin_][classification+'_id']:
                                            halo_idx = list(info_['profiles'][bin_]['all_id']).index(halo)
                                            bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                            bin_vals_.append(bin_val_)
                                elif grids[plot_kwargs['cbar']['calculate_on']]['type'] == 'bins':
                                    for bin_ in profiles_info[code]['profiles']['bins']:
                                        for halo in info['profiles'][bin_][classification+'_id']:
                                            halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                            bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                            bin_vals_.append(bin_val_)
                                elif grids[plot_kwargs['cbar']['calculate_on']]['type'] == 'props':
                                    for halo in info['profiles'][bin_][classification+'_id']:
                                        halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                        bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                        bin_vals_.append(bin_val_)
                                elif grids[plot_kwargs['cbar']['calculate_on']]['type'] == 'classifications':
                                    for classification_ in classifications:
                                        for halo in info['profiles'][bin_][classification_+'_id']:
                                            halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                            bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                            bin_vals_.append(bin_val_)


                            elif plot_kwargs['cbar']['calculate_on'] == 'all':
                                ## run over all profile binning values in both row and column
                                if grids['row']['type'] == 'sims':
                                    for code, info in profiles_info.items():
                                        if grids['col']['type'] == 'bins':
                                            for bin_ in profiles_info[code]['profiles']['bins']:
                                                for halo in info['profiles'][bin_][classification+'_id']:
                                                    halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                    bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                    bin_vals_.append(bin_val_)
                                        elif grids['col']['type'] == 'props':
                                            for halo in info['profiles'][bin_][classification+'_id']:
                                                halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)
                                        elif grids['col']['type'] == 'classifications':
                                            for classification_ in classifications:
                                                for halo in info['profiles'][bin_][classification_+'_id']:
                                                    halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                    bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                    bin_vals_.append(bin_val_)

                                elif grids['row']['type'] == 'bins':
                                    for bin_ in profiles_info[code]['profiles']['bins']:
                                        if grids['col']['type'] == 'sims':
                                            for code, info in profiles_info.items():
                                                for halo in info['profiles'][bin_][classification+'_id']:
                                                    halo_idx = list(info_['profiles'][bin_]['all_id']).index(halo)
                                                    bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                    bin_vals_.append(bin_val_)
                                        elif grids['col']['type'] == 'props':
                                            for halo in info['profiles'][bin_][classification+'_id']:
                                                halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)
                                        elif grids['col']['type'] == 'classifications':
                                            for classification_ in classifications:
                                                for halo in info['profiles'][bin_][classification_+'_id']:
                                                    halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                    bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                    bin_vals_.append(bin_val_)

                                elif grids['row']['type'] == 'props':
                                    if grids['col']['type'] == 'sims':
                                        for code, info in profiles_info.items():
                                            for halo in info['profiles'][bin_][classification+'_id']:
                                                halo_idx = list(info_['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)
                                    elif grids['col']['type'] == 'bins':
                                        for bin_ in profiles_info[code]['profiles']['bins']:
                                            for halo in info['profiles'][bin_][classification+'_id']:
                                                halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)
                                    elif grids['col']['type'] == 'classifications':
                                        for classification_ in classifications:
                                            for halo in info['profiles'][bin_][classification_+'_id']:
                                                halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)

                                elif grids['row']['type'] == 'classifications':
                                    if grids['col']['type'] == 'sims':
                                        for code, info in profiles_info.items():
                                            for halo in info['profiles'][bin_][classification+'_id']:
                                                halo_idx = list(info_['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)
                                    elif grids['col']['type'] == 'bins':
                                        for bin_ in profiles_info[code]['profiles']['bins']:
                                            for halo in info['profiles'][bin_][classification+'_id']:
                                                halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                                bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                                bin_vals_.append(bin_val_)
                                    elif grids['col']['type'] == 'props':
                                        for halo in info['profiles'][bin_][classification+'_id']:
                                            halo_idx = list(info['profiles'][bin_]['all_id']).index(halo)
                                            bin_val_ = info['profiles'][bin_]['halo_info'][plot_kwargs['binning']['prop'][0]][plot_kwargs['binning']['prop'][1]][halo_idx]
                                            bin_vals_.append(bin_val_)

                        bin_color_limits = [0.95*min(bin_vals_), 1.05*max(bin_vals_)]
                    else:
                        raise Exception("bin_color_method must be 'default', 'calculate', or 'provide'")

                    if plot_kwargs['cbar']['log_cbar']:
                        bin_color_limits = np.log10(bin_color_limits)

                    if plot_kwargs['display']['verbose']:
                        print('Bin color limits:', bin_color_limits)

                    sm = plt.cm.ScalarMappable(cmap=plot_kwargs['cbar']['cm'], 
                                               norm=mpl.colors.Normalize(vmin=bin_color_limits[0], vmax=bin_color_limits[1]))
                
                
                
                sm._A = []
                
                
                
                if plot_kwargs['cbar']['placement'].lower() == 'ind':
                    sm_arr[row,col] = sm
                elif plot_kwargs['cbar']['placement'].lower() == 'row':
                    sm_arr[row] = sm
                elif plot_kwargs['cbar']['placement'].lower() == 'col':
                    sm_arr[col] = sm
                elif plot_kwargs['cbar']['placement'].lower() == 'fig':
                    sm_arr = sm

#                 ## cbar_placement can be 'ind', 'row', 'col', or 'all', 'fig'
# #                 if not fig_cbar:

# #                     if plot_kwargs['cbar']['placement'].lower() == 'fig':
# #                         fig_cbar = True
#                 if ((plot_kwargs['cbar']['placement'].lower() == 'ind') or
#                     (plot_kwargs['cbar']['placement'].lower() == 'row' and col == ncols-1) or
#                     (plot_kwargs['cbar']['placement'].lower() == 'col' and row == 0)):
# #                     div = make_axes_locatable(main_ax)
# #                     div = make_axes_locatable(ax[row,col])
    
#                     if plot_kwargs['cbar']['placement'].lower() == 'ind':
# #                         div = make_axes_locatable(main_ax)
#                         cax, ckw = make_axes(main_ax, location=plot_kwargs['cbar']['location'],
#                                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
#                                              pad=plot_kwargs['cbar']['pad'])
#                     elif plot_kwargs['cbar']['placement'].lower() == 'row':
# #                         div = make_axes_locatable(ax[row,:].ravel())
#                         cax, ckw = make_axes(ax[row,:].ravel().tolist(), location=plot_kwargs['cbar']['location'],
#                                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
#                                              pad=plot_kwargs['cbar']['pad'])
#                     elif plot_kwargs['cbar']['placement'].lower() == 'col':
# #                         div = make_axes_locatable(ax[:,col].ravel())
#                         cax, ckw = make_axes(ax[:,col].ravel().tolist(), location=plot_kwargs['cbar']['location'],
#                                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
#                                              pad=plot_kwargs['cbar']['pad'])
    
# #                     cax = div.append_axes(position=plot_kwargs['cbar']['location'], 
# #                                           size=plot_kwargs['cbar']['width'], pad=plot_kwargs['cbar']['pad'])
                
# #                     cbar = fig.colorbar(sm, cax=cax, 
# #                                         orientation=plot_kwargs['cbar']['orientation'], spacing=plot_kwargs['cbar']['spacing'],
# #                                         drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
# #                                         fraction=plot_kwargs['cbar']['fraction'],)
        
#                     cbar = fig.colorbar(sm, cax=cax, **ckw,
#                                         spacing=plot_kwargs['cbar']['spacing'], drawedges=plot_kwargs['cbar']['drawedges'])

#                     if plot_kwargs['cbar']['label_position'].lower() == 'x':
#                         cbar.ax.set_xlabel(cbarlabel, labelpad=8, fontsize=plot_kwargs['cbar']['fontsize'])
#                     elif plot_kwargs['cbar']['label_position'].lower() == 'y':
#                         cbar.set_label(label=cbarlabel, size=plot_kwargs['cbar']['fontsize'])
#                     cbar.ax.tick_params(labelsize=plot_kwargs['cbar']['fontsize'])
                    
#                     cax.xaxis.set_ticks_position(plot_kwargs['cbar']['location'])
#                     cax.xaxis.set_label_position(plot_kwargs['cbar']['location'])
                        
#                         cbar = fig.colorbar(sm, ax=main_ax, location=plot_kwargs['cbar']['location'],
#                                             pad=plot_kwargs['cbar']['pad'], spacing=plot_kwargs['cbar']['spacing'],
#                                             drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
#                                             fraction=plot_kwargs['cbar']['fraction'], use_gridspec=plot_kwargs['cbar']['use_gridspec'])
#                     elif (plot_kwargs['cbar']['placement'].lower() == 'row' and col == ncols-1):
#                         cbar = fig.colorbar(sm, ax=main_ax, location=plot_kwargs['cbar']['location'],
#                                             pad=plot_kwargs['cbar']['pad'], spacing=plot_kwargs['cbar']['spacing'],
#                                             drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
#                                             fraction=plot_kwargs['cbar']['fraction'], use_gridspec=plot_kwargs['cbar']['use_gridspec'])
#                     elif (plot_kwargs['cbar']['placement'].lower() == 'col' and row == 0):
#                         cbar = fig.colorbar(sm, ax=main_ax, location=plot_kwargs['cbar']['location'],
#                                             pad=plot_kwargs['cbar']['pad'], spacing=plot_kwargs['cbar']['spacing'],
#                                             drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
#                                             fraction=plot_kwargs['cbar']['fraction'], use_gridspec=plot_kwargs['cbar']['use_gridspec'])
#                     elif plot_kwargs['cbar']['placement'].lower() == 'all':
# #                         and row == 0 and col == ncols-1
#                         cbar = fig.colorbar(sm, ax=main_ax, location=plot_kwargs['cbar']['location'],
#                                             pad=plot_kwargs['cbar']['pad'], spacing=plot_kwargs['cbar']['spacing'],
#                                             drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
#                                             fraction=plot_kwargs['cbar']['fraction'], use_gridspec=plot_kwargs['cbar']['use_gridspec'])

#     #                     if sim_plot_kwargs['residuals']['plot']:
#     #                         div = make_axes_locatable(residual_ax)
#     #                     else:
#     # #                         div = make_axes_locatable(ax[row,col])
#     #                         div = make_axes_locatable(main_ax)
#                         if plot_kwargs['cbar']['placement'] == 'col':
#                             div = make_axes_locatable(main_ax)
#                         else:
#                             div = make_axes_locatable(ax[row,col])
#                         if plot_kwargs['cbar']['placement'] == 'col':
#                             cax_ = div.append_axes('top', size=plot_kwargs['cbar']['width'], pad=plot_kwargs['cbar']['pad']) #'5%'
#                             orientation = 'horizontal'
#     #                         label_loc = 'top'
#     #                         label_ha = 'center'
#     #                         label_va = 'top'
#                         else:
#                             cax_ = div.append_axes('right', size=plot_kwargs['cbar']['width'], pad=plot_kwargs['cbar']['pad'])
#                             orientation = 'vertical'
#     #                         label_loc = 'right'
#     #                         label_ha = 'right'
#     #                         label_va = 'center'
#                         cbar = fig.colorbar(sm, cax=cax_, orientation=orientation)
#                         if plot_kwargs['cbar']['label_position'].lower() == 'x':
#                             cbar.ax.set_xlabel(cbarlabel, labelpad=8, fontsize=plot_kwargs['cbar']['fontsize'])#, loc=label_loc)
#     #                                            ha=label_ha, va=label_va)
#                         elif plot_kwargs['cbar']['label_position'].lower() == 'y':
#                             cbar.set_label(label=cbarlabel, size=plot_kwargs['cbar']['fontsize'])#, loc=label_loc)
#     #                                        ha=label_ha, va=label_va)
#                         cbar.ax.tick_params(labelsize=plot_kwargs['cbar']['fontsize'])
#                         if plot_kwargs['cbar']['placement'] == 'col':
#                             cax_.xaxis.set_ticks_position("top")
#                             cax_.xaxis.set_label_position("top")
            else:
                sm=None
            
            
            
#             if plot_kwargs['display']['verbose']:
#                 print()
#                 print('current_info:')
#                 print(current_info.keys())
#                 print('code:', current_info['code'])
#                 print('bin_:', current_info['bin_'])
#                 print('prop:', current_info['prop'])
#                 print('classification:', current_info['classification'])
#                 print('element:', current_info['element'])
#                 print('element_num:', current_info['element_num'])
#                 print('element_denom:', current_info['element_denom'])
#                 print()
            


            ## Observations
            if obs_plot_kwargs['plot']:
                xobs.plot_observations(current_info['prop'], current_info['info'], current_info['bin_'], type_=type_, axes=main_ax,
                                       element=current_info['element'], 
                                       element_num=current_info['element_num'], element_denom=current_info['element_denom'], 
                                       bin_low=bin_edges[current_info['bin_']], bin_high=bin_edges[current_info['bin_']+1], sm=sm,
                                       plot_kwargs=plot_kwargs, obs_plot_kwargs=obs_plot_kwargs, sim_plot_kwargs=sim_plot_kwargs)
                
            
            if sim_plot_kwargs['plot']:
                num_loops = 0
                num_per_loop = []
                for key, grid_type in grids.items():
                    if 'overlay' in key:
                        num_loops += 1
                        num_per_loop.append(grid_type['length'])
                        
                overlay_args = [current_info]#, current_axis_info]
                overlay_kwargs = {
                    'grids':grids,
                    'all_grid_types':all_grid_types,
                    'bin_info':bin_info,
                    'profiles_info':profiles_info,
                    'props':props,
                    'type_':type_,
                    'elements':elements,
                    'elements_num':elements_num,
                    'elements_denom':elements_denom,
                    'plot_kwargs':plot_kwargs,
                    'sim_plot_kwargs':sim_plot_kwargs,
                    'current_axis_info':current_axis_info,
                }
                
#                 print(current_axis_info)
                plot_overlay_args = [current_info]#, current_axis_info]
                plot_overlay_kwargs = {
                    'current_axis_info':current_axis_info,
                    'ax':main_ax,
                    'sm':sm,
                    'profiles_info':profiles_info,
#                     'props':props,
                    'type_':type_,
#                     'elements':elements,
#                     'elements_num':elements_num,
#                     'elements_denom':elements_denom,
                    'plot_kwargs':plot_kwargs,
                    'sim_plot_kwargs':sim_plot_kwargs,
                }
                
                ## Plot property
                recursive_loop(num_loops,
                               num_per_loop=num_per_loop,
                               func=update_info_and_create_label, trigger_func=plot_overlay, 
                               args=overlay_args, kwargs=overlay_kwargs, 
                               trigger_args=plot_overlay_args, trigger_kwargs=plot_overlay_kwargs)
                
                
                ## Actually plot residuals
                if sim_plot_kwargs['residuals']['plot']:
                    all_main_lines, all_main_labels = main_ax.get_legend_handles_labels()
                    main_lines, main_labels = [], []
                    for main_line, main_label in zip(all_main_lines, all_main_labels):
#                         if 'avg!' in main_label and 'sim' in main_label:
                        if np.all([key in main_label for key in sim_plot_kwargs['residuals']['need_in_label']]):
                            main_labels.append(main_label)
                            main_lines.append(main_line)
#                     print('main_lines:')
#                     print(main_lines)
#                     print()
#                     print('main_labels:')
#                     print(main_labels)
#                     print()
                    for ll in range(len(sim_plot_kwargs['residuals']['compare'])):
#                         print()
#                         print('sim_plot_kwargs[residuals][compare]:', sim_plot_kwargs['residuals']['compare'])
#                         print()
                        plotted_sim_idxs = sim_plot_kwargs['residuals']['compare'][ll]
#                         print()
#                         print('plotted_sim_idxs:', plotted_sim_idxs)
#                         print()
                        
                        ## Get labels and lines for sim_profiles from main_ax
                        idx_num = plotted_sim_idxs[0]
                        idx_denom = plotted_sim_idxs[1]
                        
                        try:
                            x_num, y_num = main_lines[idx_num].get_data()
                            x_denom, y_denom = main_lines[idx_denom].get_data()
                        except:
                            print('Not plotting residuals; likely one of the datasets does not exist')
                            continue
                        
                        assert np.array_equal(x_num, x_denom, equal_nan=True) , '%s and %s simulation profiles have different x-axes. Cannot calculate residuals.' % (label_num, label_denom)
                        
                        ## Calculate residual according to provided method
                        residual = sim_plot_kwargs['residuals']['method'](y_num, y_denom)
                        
                        ## determine label
                        if sim_plot_kwargs['residuals']['provide_labels']:
                            residual_label = sim_plot_kwargs['residuals']['labels'][ll]
                        else:
                            label_num = main_labels[idx_num]
                            label_denom = main_labels[idx_denom]
#                             print('label_num:', label_num)
#                             print('label_denom:', label_denom)
#                             print()
                            
                            if '##' in label_num:
                                label_num = label_num.split('##')[-1]
                            if '##' in label_denom:
                                label_denom = label_denom.split('##')[-1]
#                             print('label_num:', label_num)
#                             print('label_denom:', label_denom)
#                             print()
                            
                            residual_label = 'avg!#sim##' + label_num + ' / ' + label_denom
#                             print('label_num:', label_num)
#                             print('label_denom:', label_denom)
#                         print()
#                         print('residual_label:', residual_label)
#                         print()
                        
                        residual_ax.plot(gs_xaxis, residual, 
                                     color=sim_plot_kwargs['residuals']['color_list'][ll],
                                     ls=sim_plot_kwargs['residuals']['ls'][ll], 
                                     lw=sim_plot_kwargs['residuals']['lw'][ll], 
                                     alpha=sim_plot_kwargs['residuals']['alpha'], 
                                     zorder=sim_plot_kwargs['residuals']['zorder'], 
                                     label=residual_label)
                        
                    

            
            if plot_kwargs['display']['verbose']:
                print()
            


# #             for axis_ in [main_ax, residual_ax]:
# #                 plot_extra_lines(axis=axis_, xunits=current_info['info']['XAXIS'].units, plot_kwargs=plot_kwargs)
#             plot_extra_lines(axis=main_ax, xunits=current_info['info']['XAXIS'].units, plot_kwargs=plot_kwargs)
#             if sim_plot_kwargs['residuals']['plot']:
#                 plot_extra_lines(axis=residual_ax, xunits=current_info['info']['XAXIS'].units, plot_kwargs=plot_kwargs)
                
                
#             plot_extra_lines_new(axis=main_ax, plot_kwargs=plot_kwargs['extra_lines'])
#             if sim_plot_kwargs['residuals']['plot']:
#                 plot_extra_lines_new(axis=residual_ax, plot_kwargs=plot_kwargs['residuals_extra_lines'])


            
            ## Get rid of residual subplot only for current subplot if there are no lines in residual subplot
            add_main_label = False
            if sim_plot_kwargs['residuals']['plot']:
                residual_has_lines = True
                residual_lines = residual_ax.get_lines()
                if len(residual_lines) == 0:
                    residual_has_lines = False
                    
                    ## Copy contents of main_ax to add to a new subplot that takes up the whole space of the subplot
# #                     main_lines_copy, main_labels_copy = copy.deepcopy(main_ax.get_legend_handles_labels())
#                     main_lines_copy, main_labels_copy = main_ax.get_legend_handles_labels()
    
# #                     gs_outer[row,col].remove()
# #                     main_ax.remove()
# #                     residual_ax.remove()
#                     main_ax.cla()
#                     residual_ax.cla()
#                     main_ax.clear()
#                     residual_ax.clear()
# #                     main_ax.set_axis_off()
# #                     residual_ax.set_axis_off()
#                     main_ax.remove()
#                     residual_ax.remove()
                    
# #                     gs_inner.remove()
                    
# #                     main_x_copy, main_y_copy = [], []
# #                     for main_line in main_lines_copy:
# #                         main_x, main_y = main_lines.get_data()
# #                         main_x_copy.append(main_x)
# #                         main_y_copy.append(main_y)
                    
#                     main_ax = ax[row,col]
#                     residual_ax = None

#                     if plot_kwargs['figure']['sharex'] != False and plot_kwargs['figure']['sharex'] != 'none' and row == nrows-1:
#                         main_ax.tick_params(labelbottom=True)


#                     if (plot_kwargs['figure']['sharex'] != False and 
#                         plot_kwargs['figure']['sharex'] != 'none' and 
#                         plot_kwargs['figure']['sharex'] != 'row'):                        
#                         if row == nrows-1:
#                             main_ax.tick_params(labelbottom=True)
#                         else:
#                             main_ax.tick_params(labelbottom=False)
                            
#                     for main_line in main_lines_copy:
#                         main_ax.add_line(main_line)
                    
    
                    ax[row,col].remove()
#                     gs_inner.set_visible(False)
#                     residual_ax.axis('off')
#                     residual_ax.set_visible(False)
#                     residual_ax.cla()
#                     residual_ax.clear()
                    residual_ax.remove()
#                     gs_inner[1].remove()

                    if row == nrows-1:
                        add_main_label = True
#                         print()
#                         print('nrows:',nrows)
#                         print('last row:',row)
#                         print()
#                         print('HIIIIIII')
                        main_ax.axis('on')
#                         main_ax.tick_params(labelbottom=True)
#                         ax[row,col].tick_params(labelbottom=True)
            

                
            if plot_kwargs['extra_lines']['plot']:
                extra_line_idx = row if plot_kwargs['extra_lines']['plot_along'].lower() == 'row' else col
                plot_extra_lines_new(labels=plot_kwargs['extra_lines']['labels'][extra_line_idx],
                                     orientations=plot_kwargs['extra_lines']['orientations'][extra_line_idx], 
                                     values=plot_kwargs['extra_lines']['values'][extra_line_idx], 
                                     ls=plot_kwargs['extra_lines']['ls'][extra_line_idx],
                                     lw=plot_kwargs['extra_lines']['lw'][extra_line_idx],
                                     colors=plot_kwargs['extra_lines']['colors'][extra_line_idx],
                                     zorders=plot_kwargs['extra_lines']['zorders'][extra_line_idx], 
                                     axis=main_ax)
            
            if plot_kwargs['residuals_extra_lines']['plot'] and sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                residuals_extra_line_idx = row if plot_kwargs['residuals_extra_lines']['plot_along'].lower() == 'row' else col
                plot_extra_lines_new(labels=plot_kwargs['residuals_extra_lines']['labels'][residuals_extra_line_idx],
                                     orientations=plot_kwargs['residuals_extra_lines']['orientations'][residuals_extra_line_idx], 
                                     values=plot_kwargs['residuals_extra_lines']['values'][residuals_extra_line_idx], 
                                     ls=plot_kwargs['residuals_extra_lines']['ls'][residuals_extra_line_idx],
                                     lw=plot_kwargs['residuals_extra_lines']['lw'][residuals_extra_line_idx],
                                     colors=plot_kwargs['residuals_extra_lines']['colors'][residuals_extra_line_idx],
                                     zorders=plot_kwargs['residuals_extra_lines']['zorders'][residuals_extra_line_idx], 
                                     axis=residual_ax)


            if row == 0 and plot_kwargs['col_titles']['show']:
#                 ax[0,col].set_title(col_label, wrap=plot_kwargs['col_titles']['wrap'], 
#                                     fontsize=plot_kwargs['col_titles']['fontsize'])
                main_ax.set_title(col_label, wrap=plot_kwargs['col_titles']['wrap'], 
                                    fontsize=plot_kwargs['col_titles']['fontsize'])

            if row == nrows-1:
                if add_main_label:
                    main_ax.set_xlabel(plot_kwargs['xaxis']['label'], 
                                       fontsize=plot_kwargs['xaxis']['label_fontsize'], 
                                       labelpad=plot_kwargs['xaxis']['labelpad'], zorder=999)
                else:
                    ax[row,col].set_xlabel(plot_kwargs['xaxis']['label'], 
                                       fontsize=plot_kwargs['xaxis']['label_fontsize'], 
                                       labelpad=plot_kwargs['xaxis']['labelpad'])
#                 ax[row,col].set_xlabel(plot_kwargs['xaxis']['label'], fontsize=plot_kwargs['xaxis']['label_fontsize'])
            
            
            if col == 0:
                if plot_kwargs['yaxis']['use_alt_yaxis_labels']:
                    ylabel_ = plot_kwargs['yaxis']['alt_yaxis_labels'][row]
                else:
                    ylabel_ = plot_kwargs['prop_specific_names'][row]
                    
                if plot_kwargs['yaxis']['show_extra_info']:
                    if type_ == 'metallicity':
                        ylabel_ += ' [%s]' % element
                    elif type_ == 'metallicity ratio':
                        ylabel_ += ' [%s/%s]' % (element_num, element_denom)
                    
                if plot_kwargs['yaxis']['is_log']:
                    ylabel_ = r'$\log($%s$)$' % ylabel_
                
                if plot_kwargs['grid']['show_prop_labels']:
                    ylabel_ = '%s \n \n %s' % (row_label, ylabel_)
                    
#                 ax[row,col].set_ylabel(ylabel_, fontsize=plot_kwargs['yaxis']['label_fontsize'])
                main_ax.set_ylabel(ylabel_, fontsize=plot_kwargs['yaxis']['label_fontsize'])
    
                if sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                    residual_ax.set_ylabel(plot_kwargs['yaxis']['residuals']['label'],
                                       fontsize=plot_kwargs['yaxis']['residuals']['label_fontsize'])


            if plot_kwargs['xaxis']['plot_log']:
#                 ax[row,col].set_xscale('log')
                main_ax.set_xscale('log')
                if sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                    residual_ax.set_xscale('log')

            if plot_kwargs['yaxis']['plot_log']:
#                 ax[row,col].set_yscale('log')
                main_ax.set_yscale('log')
            if sim_plot_kwargs['residuals']['plot'] and residual_has_lines and plot_kwargs['yaxis']['residuals']['plot_log']:
                residual_ax.set_yscale('log')
            
            if plot_kwargs['xaxis']['set_lims']:
#                 ax[row,col].set_xlim(plot_kwargs['xaxis']['lims'][0], plot_kwargs['xaxis']['lims'][1])
                main_ax.set_xlim(plot_kwargs['xaxis']['lims'][0], plot_kwargs['xaxis']['lims'][1])
                if sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                    residual_ax.set_xlim(plot_kwargs['xaxis']['lims'][0], plot_kwargs['xaxis']['lims'][1])
        
            if plot_kwargs['yaxis']['set_lims']:
#                 ax[row,col].set_ylim(plot_kwargs['yaxis']['lo_lims'][row], plot_kwargs['yaxis']['hi_lims'][row])
                main_ax.set_ylim(plot_kwargs['yaxis']['lo_lims'][row], plot_kwargs['yaxis']['hi_lims'][row])
            
            if sim_plot_kwargs['residuals']['plot'] and residual_has_lines and plot_kwargs['yaxis']['residuals']['set_lims']:
                residual_ax.set_ylim(plot_kwargs['yaxis']['residuals']['lo_lims'][row], plot_kwargs['yaxis']['residuals']['hi_lims'][row])

#             ax[row,col].minorticks_on()
            main_ax.minorticks_on()
            if sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                residual_ax.minorticks_on()
#                 ax[row,col].minorticks_off()
            if plot_kwargs['figure']['show_grid']:
#                 ax[row,col].grid()
                main_ax.grid()
                if sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                    residual_ax.grid()

#             ax[row,col].tick_params(axis='x', **plot_kwargs['ticks'])
#             ax[row,col].tick_params(axis='y', **plot_kwargs['ticks'])
            main_ax.tick_params(axis='x', **plot_kwargs['ticks'])
            main_ax.tick_params(axis='y', **plot_kwargs['ticks'])
            if sim_plot_kwargs['residuals']['plot'] and residual_has_lines:
                residual_ax.tick_params(axis='x', **plot_kwargs['ticks'])
                residual_ax.tick_params(axis='y', **plot_kwargs['ticks'])
                
                ax[row,col].tick_params(axis='x', colors='none')
                ax[row,col].tick_params(axis='y', colors='none')
                
                


                
    ## Add extra text to plots
    for text_info in plot_kwargs['extra_text'].values():
#         print(text_info)
        if sim_plot_kwargs['residuals']['plot']:
#             print('HI')
#             print(text_info['outer_axis'][0])
#             print(text_info['outer_axis'][1])
#             print(text_info['inner_axis'])
            text_axis = gs_axis_list[text_info['outer_axis'][0]][text_info['outer_axis'][1]][text_info['inner_axis']]
        else:
#             print('BYE')
#             text_axis = ax[text_info['outer_axis']]
            text_axis = ax[text_info['outer_axis'][0]][text_info['outer_axis'][1]]
            
#         print(ax)
#         print()
#         print(text_axis)

        if text_info['coords'].lower() == 'data':
            text_axis.text(**text_info['dict'], transform=text_axis.transData)
        elif text_info['coords'].lower() == 'axes':
            text_axis.text(**text_info['dict'], transform=text_axis.transAxes)
        else:
            text_axis.text(**text_info['dict'])
            
            
    ## Add cbar to whole figure
#     if fig_cbar:

    if plot_kwargs['cbar']['plot']:
        if plot_kwargs['cbar']['placement'].lower() == 'ind':
            for row in range(nrows):
                for col in range(ncols):
                    if sim_plot_kwargs['residuals']['plot']:
    #                     cbar_ax = gs_axis_dict['main_ax'][row,col]
                        cbar_ax = gs_axis_dict['outer_ax'][row,col]
                    else:
                        cbar_ax = ax[row,col]
                    cax, ckw = make_axes(cbar_ax, location=plot_kwargs['cbar']['location'],
                                         shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
                                         pad=plot_kwargs['cbar']['pad'])

                    sm = sm_arr[row,col]
                    cbar = fig.colorbar(sm, cax=cax, **ckw,
                                        spacing=plot_kwargs['cbar']['spacing'], drawedges=plot_kwargs['cbar']['drawedges'])

                    if plot_kwargs['cbar']['label_position'].lower() == 'x':
                        cbar.ax.set_xlabel(cbarlabel, labelpad=8, fontsize=plot_kwargs['cbar']['fontsize'])
                    elif plot_kwargs['cbar']['label_position'].lower() == 'y':
                        cbar.set_label(label=cbarlabel, size=plot_kwargs['cbar']['fontsize'])
                    cbar.ax.tick_params(labelsize=plot_kwargs['cbar']['fontsize'])

                    cax.xaxis.set_ticks_position(plot_kwargs['cbar']['location'])
                    cax.xaxis.set_label_position(plot_kwargs['cbar']['location'])

        elif plot_kwargs['cbar']['placement'].lower() == 'row':
            for row in range(nrows):
                if sim_plot_kwargs['residuals']['plot']:
                    cbar_ax = gs_axis_dict['main_ax'][row,:]
                else:
                    cbar_ax = ax[row,:]
                cax, ckw = make_axes(cbar_ax, location=plot_kwargs['cbar']['location'],
                                     shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
                                     pad=plot_kwargs['cbar']['pad'])

                sm = sm_arr[row]
                cbar = fig.colorbar(sm, cax=cax, **ckw,
                                    spacing=plot_kwargs['cbar']['spacing'], drawedges=plot_kwargs['cbar']['drawedges'])

                if plot_kwargs['cbar']['label_position'].lower() == 'x':
                    cbar.ax.set_xlabel(cbarlabel, labelpad=8, fontsize=plot_kwargs['cbar']['fontsize'])
                elif plot_kwargs['cbar']['label_position'].lower() == 'y':
                    cbar.set_label(label=cbarlabel, size=plot_kwargs['cbar']['fontsize'])
                cbar.ax.tick_params(labelsize=plot_kwargs['cbar']['fontsize'])

                cax.xaxis.set_ticks_position(plot_kwargs['cbar']['location'])
                cax.xaxis.set_label_position(plot_kwargs['cbar']['location'])
        elif plot_kwargs['cbar']['placement'].lower() == 'col':
            for col in range(ncols):
                if sim_plot_kwargs['residuals']['plot']:
                    cbar_ax = gs_axis_dict['main_ax'][:,col]
                else:
                    cbar_ax = ax[:,col]
                cax, ckw = make_axes(cbar_ax, location=plot_kwargs['cbar']['location'],
                                     shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
                                     pad=plot_kwargs['cbar']['pad'])

                sm = sm_arr[col]
                cbar = fig.colorbar(sm, cax=cax, **ckw,
                                    spacing=plot_kwargs['cbar']['spacing'], drawedges=plot_kwargs['cbar']['drawedges'])

                if plot_kwargs['cbar']['label_position'].lower() == 'x':
                    cbar.ax.set_xlabel(cbarlabel, labelpad=8, fontsize=plot_kwargs['cbar']['fontsize'])
                elif plot_kwargs['cbar']['label_position'].lower() == 'y':
                    cbar.set_label(label=cbarlabel, size=plot_kwargs['cbar']['fontsize'])
                cbar.ax.tick_params(labelsize=plot_kwargs['cbar']['fontsize'])

                cax.xaxis.set_ticks_position(plot_kwargs['cbar']['location'])
                cax.xaxis.set_label_position(plot_kwargs['cbar']['location'])

        elif plot_kwargs['cbar']['placement'].lower() == 'fig':


    #         cax = add_cax_to_figure(fig, nrows=nrows, ncols=ncols, location=plot_kwargs['cbar']['location'])
    #         cax = add_cax_to_figure(fig, nrows=ax.shape[0], ncols=ax.shape[1], location=plot_kwargs['cbar']['location'])
            cax = add_cax_to_figure(fig, fig, ax_is_fig=True, x=plot_kwargs['cbar']['x'], y=plot_kwargs['cbar']['y'],
                                    width=plot_kwargs['cbar']['width'], height=plot_kwargs['cbar']['height'],
                                    x_in_inches=plot_kwargs['cbar']['x_in_inches'], y_in_inches=plot_kwargs['cbar']['y_in_inches'],
                                    width_in_inches=plot_kwargs['cbar']['width_in_inches'],
                                    height_in_inches=plot_kwargs['cbar']['height_in_inches'])

    #         sm = sm_arr
    #         cbar = 




    #         if sim_plot_kwargs['residuals']['plot']:
    #             cbar_ax = gs_axis_dict['inner_ax']
    #         else:
    #             cbar_ax = ax

    #         cax, ckw = make_axes(cbar_ax, #location=plot_kwargs['cbar']['location'],
    #                              orientation=plot_kwargs['cbar']['orientation'],
    #                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
    #                              pad=plot_kwargs['cbar']['pad'])

            sm = sm_arr
            cbar = fig.colorbar(sm, cax=cax, #location=plot_kwargs['cbar']['location'], #**ckw, # 
                                orientation=plot_kwargs['cbar']['orientation'],
                                spacing=plot_kwargs['cbar']['spacing'], drawedges=plot_kwargs['cbar']['drawedges'],
    #                             ticklocation=plot_kwargs['cbar']['location']
                                shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
                                pad=plot_kwargs['cbar']['pad'],
                               )

            if plot_kwargs['cbar']['label_position'].lower() == 'x':
                cbar.ax.set_xlabel(cbarlabel, labelpad=8, fontsize=plot_kwargs['cbar']['fontsize'])
            elif plot_kwargs['cbar']['label_position'].lower() == 'y':
                cbar.set_label(label=cbarlabel, size=plot_kwargs['cbar']['fontsize'])
            cbar.ax.tick_params(labelsize=plot_kwargs['cbar']['fontsize'])

            if plot_kwargs['cbar']['location'].lower() in ['left', 'right']:
                cax.yaxis.set_ticks_position(plot_kwargs['cbar']['location'])
                cax.yaxis.set_label_position(plot_kwargs['cbar']['location'])
            if plot_kwargs['cbar']['location'].lower() in ['top', 'bottom']:
                cax.xaxis.set_ticks_position(plot_kwargs['cbar']['location'])
                cax.xaxis.set_label_position(plot_kwargs['cbar']['location'])
            
            
            
            
        
# #         div = make_axes_locatable(ax.flat)
# #         div = make_axes_locatable(fig)
# #         cax = fig.append_axes(position=plot_kwargs['cbar']['location'], 
# #                               size=plot_kwargs['cbar']['width'],
# #                               pad=plot_kwargs['cbar']['pad'])
# #         cax = fig.add_axes(position=plot_kwargs['cbar']['location'], 
# #                               size=plot_kwargs['cbar']['width'],
# #                               pad=plot_kwargs['cbar']['pad'])
# #         cax, ckw = make_axes([axis for axis in ax.flat], location=plot_kwargs['cbar']['location'],
# #                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
# #                              pad=plot_kwargs['cbar']['pad'])
#         if sim_plot_kwargs['residuals']['plot']:
# #             fig_ax = gs_outer#.ravel().tolist()
# #             fig_ax = gs_axis_dict['inner_ax'].ravel().tolist()
# #             cbar_ax = gs_axis_dict['outer_ax']#.ravel().tolist()
# #             cbar_ax = gs_axis_dict['main_ax']#.ravel().tolist()
#             cbar_ax = gs_axis_dict['inner_ax']
# #             print(fig_ax)
# #             cax, ckw = make_axes_gridspec(fig_ax, location=plot_kwargs['cbar']['location'],
# #                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
# #                              pad=plot_kwargs['cbar']['pad'])
# #             cax, ckw = make_axes(gs_axis_dict['all_ax'].ravel().tolist(), location=plot_kwargs['cbar']['location'], 
# # #                                           orientation=plot_kwargs['cbar']['orientation'],
# #                                  shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
# #                                  pad=plot_kwargs['cbar']['pad'])
#         else:
# #             fig_ax = ax.ravel().tolist()
#             cbar_ax = ax#.ravel().tolist()
# #             cax, ckw = make_axes(ax.ravel().tolist(), location=plot_kwargs['cbar']['location'],
# #                                  shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
# #                                  pad=plot_kwargs['cbar']['pad'])
#         cax, ckw = make_axes(cbar_ax, #location=plot_kwargs['cbar']['location'],
#                              orientation=plot_kwargs['cbar']['orientation'],
#                              shrink=plot_kwargs['cbar']['shrink'], fraction=plot_kwargs['cbar']['fraction'],
#                              pad=plot_kwargs['cbar']['pad'])
# #         orientation=plot_kwargs['cbar']['orientation'], 
    
    
#         sm = sm_arr
#         cbar = fig.colorbar(sm, cax=cax, **ckw, #location=plot_kwargs['cbar']['location'], 
# #                             orientation=plot_kwargs['cbar']['orientation'],
#                             spacing=plot_kwargs['cbar']['spacing'], drawedges=plot_kwargs['cbar']['drawedges'],
# #                             ticklocation=plot_kwargs['cbar']['location']
#                            )
        
# #         cbar = fig.colorbar(sm, cax=cax, 
# #                             orientation=plot_kwargs['cbar']['orientation'], spacing=plot_kwargs['cbar']['spacing'],
# #                             drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
# #                             fraction=plot_kwargs['cbar']['fraction'],)
        
        
        
# #         cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
# #         plt.colorbar(im, cax=cax, **kw)
        
        
# # #         if sim_plot_kwargs['residuals']['plot']:
# # #             fig_axes = np.array(gs_axis_list).flat
# # #         else:
# # #             fig_axes = 
# #         cbar = fig.colorbar(sm, ax=ax.flat, location=plot_kwargs['cbar']['location'],
# #                             pad=plot_kwargs['cbar']['pad'], spacing=plot_kwargs['cbar']['spacing'],
# #                             drawedges=plot_kwargs['cbar']['drawedges'], shrink=plot_kwargs['cbar']['shrink'],
# #                             fraction=plot_kwargs['cbar']['fraction'], use_gridspec=plot_kwargs['cbar']['use_gridspec'],
# #                             fontsize=plot_kwargs['cbar']['fontsize'])

                
            
    del current_info
            

    ## Legend
    if plot_kwargs['legend']['show']:

        # Define a dictionary to store legend information    
        legend_info = {}

        # Get the legend information from all subplots
#         if sim_plot_kwargs['residuals']['plot']:
#             for ax_ in gs_outer.flatten():
#                 lines, labels = ax_.get_legend_handles_labels()
#                 for line, label in zip(lines, labels):
#                     legend_info[label] = line  # Store line information with label
#         else:
#             for ax_ in ax.flatten():
#                 lines, labels = ax_.get_legend_handles_labels()
#                 for line, label in zip(lines, labels):
#                     legend_info[label] = line  # Store line information with label
        for ax_ in np.array(fig.axes).flatten():
            lines, labels = ax_.get_legend_handles_labels()
            for line, label in zip(lines, labels):
                legend_info[label] = line  # Store line information with label
                  
                  

        # Create a new subplot above the existing ones for the legend
        plt.subplots_adjust(top=0.85)  # Adjust the top margin to make space for the legend
        plt.subplots_adjust(bottom=0.15)  # Adjust the top margin to make space for the legend
        legend_ax = fig.add_subplot(111, frameon=False)
        legend_ax.get_xaxis().set_visible(False)
        legend_ax.get_yaxis().set_visible(False)

        # Add the legend with the unique labels, colors, and formats
        handles = [legend_info[label] for label in legend_info.keys()]
        labels = list(legend_info.keys())
        
        ## If only want to show extra labels
        if plot_kwargs['legend']['only_show_extra']:
            labels = []
            handles = []
        
        ## Manually add extra specified legend entries
        for label, handle in plot_kwargs['legend']['extra'].items():
            labels.append(label)
            handles.append(handle)
        
        for ll in range(len(labels)):
            label = labels[ll]
            if '!' in label:
                labels[ll] = label.split('!')[1]
        
        
        ## Combine handles for data that has, for example, both a line and a shaded region
        ## acting as its errorbar, so that the legend handle has both
        if plot_kwargs['legend']['combine_handles']:
            labels_copy = labels.copy()
            handles_copy = handles.copy()
            
            for ii in range(len(labels_copy)):
                label = labels_copy[ii]
                handle = handles_copy[ii]
                
                if label in labels_copy[:ii]+labels_copy[ii+1:] and label!=None:
                    handles[ii] = tuple(handle_ for label_, handle_ in zip(labels, handles) if label_==label)
                    handles, labels = set_list_duplicates_to_none(handles, sort_list=labels, 
                                                                  all_duplicates=False, specific_value=label,
                                                                  remove_followers=True, return_sort_list=True)
                    labels_copy[ii] = None
                    handles_copy[ii] = None
                
            labels = [label for label in labels if label!=None]
            handles = [handle for handle in handles if handle!=None]
                    

        # Order labels according to provided order
        if plot_kwargs['legend']['order']:
            label_flags = []
            # Find label flags
            for label_ in labels:
                if label_[0] == '#':
                    label_flag = label_.split('##')[0][1:]
                    if label_flag not in plot_kwargs['legend']['label_order']:
                        label_flags.append('other')
                    else:
                        label_flags.append(label_flag)
                else:
                    label_flags.append('other')

            ordered_label_flags = []  ## Provided ordered list of flags, but now with the total number of labels plus 'others'
            for flag in plot_kwargs['legend']['label_order']:
                for num in range(label_flags.count(flag)):
                    ordered_label_flags.append(flag)
            for num in range(label_flags.count('other')):
                ordered_label_flags.append('other')

            # Get list of indices of label_flags in ordered_label_flags (ie. order that labels should go in)
            ordered_idxs = []
            for i in range(len(label_flags)):
                idx = ordered_label_flags.index(label_flags[i])
                ordered_label_flags[idx] = None
                ordered_idxs.append(idx)

            ## Sort by index from ordered_idxs
            combined_idxs_labels_handles = zip(ordered_idxs, labels, handles)
            sorted_idxs_labels_handles = sorted(combined_idxs_labels_handles, key=lambda val:val[0])

            idxs = [idx for idx, label, handle in sorted_idxs_labels_handles]
            labels = [label for idx, label, handle in sorted_idxs_labels_handles]
            handles = [handle for idx, label, handle in sorted_idxs_labels_handles]
            
            
        if plot_kwargs['legend']['make_columns']:
#             if not plot_kwargs['legend']['use_ncol']:
#                 plot_kwargs['legend']['ncol'] = len(plot_kwargs['legend']['columns'].keys())
            
            legend_col_dict = {
                col_name:{
                    'labels':[], 
                    'handles':[],
                } for col_name in plot_kwargs['legend']['columns'].values()
            }
            max_col_length = 0
            for col_shortname, col_name in plot_kwargs['legend']['columns'].items():
#                 print('col_name:', col_name)
                for label, handle in zip(labels, handles):
                    if col_shortname in label:
                        legend_col_dict[col_name]['labels'].append(label)
                        legend_col_dict[col_name]['handles'].append(handle)
                        
#                 print('len(legend_col_dict[col_name][labels]):', len(legend_col_dict[col_name]['labels']))
                
                max_col_length = max(max_col_length, len(legend_col_dict[col_name]['labels']))
#                 print('max_col_length:', max_col_length)
                max_col_length = min(max_col_length, plot_kwargs['legend']['max_col_length'])
#                 print('max_col_length:', max_col_length)
                
                ## If number of labels/handles in column is greater than manually provided max_col_length,
                ## push them into next column, by creating new column with similar(?) title
                max_length_diff = len(legend_col_dict[col_name]['labels']) - plot_kwargs['legend']['max_col_length']
#                 print('max_length_diff:', max_length_diff)
#                 print()
                if max_length_diff > 0:
                    new_col_name = col_name+'...'
                    legend_col_dict[new_col_name] = {'labels':[], 'handles':[]}
#                     excess_keys = list(legend_col_dict[col_name].keys())[-1*int(max_length_diff):]
#                     excess_start_idxs = legend_col_dict[col_name]['labels']
                    for excess_idx in range(plot_kwargs['legend']['max_col_length'], len(legend_col_dict[col_name]['labels'])):
#                         print(excess_idx)
#                         excess_label = legend_col_dict[col_name]['labels'].pop(excess_idx)
#                         excess_handle = legend_col_dict[col_name]['handles'].pop(excess_idx)
                        excess_label = legend_col_dict[col_name]['labels'].pop(plot_kwargs['legend']['max_col_length'])
                        excess_handle = legend_col_dict[col_name]['handles'].pop(plot_kwargs['legend']['max_col_length'])
                        legend_col_dict[new_col_name]['labels'].append(excess_label)
                        legend_col_dict[new_col_name]['handles'].append(excess_handle)
#                     for excess_key in excess_keys:
#                         legend_col_dict['also ' + col_name]['labels'].append()

            if not plot_kwargs['legend']['use_ncol']:
#                 plot_kwargs['legend']['ncol'] = len(plot_kwargs['legend']['columns'].keys())
                plot_kwargs['legend']['ncol'] = len(legend_col_dict.keys())
                
            labels = []
            handles = []
            for col_name, info in legend_col_dict.items():
#                 print()
#                 print('col_name:', col_name)
                length_diff = max_col_length - len(info['labels'])
                
                if plot_kwargs['legend']['use_ncol']:
                    labels += [col_name] + [label for label in info['labels']] + ['']
                    handles += [mlines.Line2D([],[],color='none',label=col_name)] + \
                                [handle for handle in info['handles']] + \
                                [mlines.Line2D([],[],color='none')]
                else:
                    labels += [col_name] + [label for label in info['labels']] + [''] * length_diff
                    handles += [mlines.Line2D([],[],color='none',label=col_name)] + \
                                [handle for handle in info['handles']] + \
                                [mlines.Line2D([],[],color='none')] * length_diff
                    
                            
        for jj in range(len(labels)):
            if len(labels[jj]) > 0:
                if labels[jj][0] == '#':
                    labels[jj] = labels[jj].split('##')[1]
                    

        if len(labels) > 0:
            
            legend = plt.legend(handles=handles, labels=labels, loc='upper center', ncol=plot_kwargs['legend']['ncol'], 
                                bbox_to_anchor=(plot_kwargs['legend']['x'], plot_kwargs['legend']['y']), 
                                **plot_kwargs['legend']['kwargs'])
            legend_ax.add_artist(legend)
            


    # Save and show plot
    if plot_kwargs['display']['save']:
        print('Saving...')
        plt.savefig(os.path.join('plots', '%s.%s' % (plot_kwargs['display']['save_name'], plot_kwargs['display']['save_format'])), 
                    format=plot_kwargs['display']['save_format'], bbox_inches='tight', 
                    facecolor=plot_kwargs['display']['facecolor'], edgecolor=plot_kwargs['display']['edgecolor'],
                    dpi=plot_kwargs['display']['dpi'], transparent=plot_kwargs['display']['transparent'])
        print('done')

    if plot_kwargs['display']['show']:
        plt.show()