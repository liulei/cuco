#ifndef	ALLVARS_H
#include	"allvars.h"
#endif

void	read_parameter_file(void);
void	set_units(void);
void	init(void);
void	check_omega(void);
void	read_ic(char * fname);
void	set_softenings(void);
void	run(void);
void	begrun(void);
void	long_range_init(void);
void	init_drift_table(void);
double	drift_integ(double a, void * param);
double	gravkick_integ(double a, void * param);
void	pm_init_periodic(void);
void	move_particles(int time0, int time1);
void	do_box_wrapping(void);
void	find_next_sync_point_and_drift(void);
int		find_next_outputtime(int ti_curr);
void	every_timestep_stuff(void);
void	build_linklist(void);
void	compute_accelerations(void);
void	advance_and_find_timesteps(void);
double	get_drift_factor(int time0, int time1);
double	get_gravkick_factor(int time0, int time1);
void	savepositions(int num);
int		read_outputlist(char * fname);
void	long_range_force(void);
void	pm_init_periodic_allocate(void);
void	pmforce_periodic(void);
void	pm_init_periodic_free(void);
int		get_timestep(int p);

#ifdef	LINKLIST
void	gravity_linklist(void);
void	linklist_init(void);
void	force_evaluate_shortrange(int target);
#endif

#ifdef	TREE
void	gravity_tree(void);
void	force_treeallocate(int maxnodes, int maxpart);
void	force_treebuild(int npart);
void	force_treeevaluate_shortrange(int target);
void	force_update_node_recursive(int no, int sib, int father);
#endif
