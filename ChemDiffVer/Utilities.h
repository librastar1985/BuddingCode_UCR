#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "SystemStructures.h"
#include <list>
#include <cmath>

class Utilities {
    std::vector<bool> boundary_node;
    std::vector<int> nndata;

    public:
    Utilities(CoordInfoVecs& coordInfoVecs, GeneralParams& generalParams);

	std::vector<bool> DomainBd (CoordInfoVecs& coordInfoVecs);
	std::vector<int> Number_of_Neighbor(CoordInfoVecs& coordInfoVecs);

    int growth_host_vecs(int iedge,
    	GeneralParams& generalParams,
	HostSetInfoVecs& hostSetInfoVecs,
	CoordInfoVecs& coordInfoVecs,
	LinearSpringInfoVecs& linearSpringInfoVecs,
	BendingTriangleInfoVecs& bendingTriangleInfoVecs,
	AreaTriangleInfoVecs& areaTriangleInfoVecs);

    int surfaceNormal_device_vecs(int inode,
        CoordInfoVecs& coordInfoVecs,
        GeneralParams& generalParams
    );

    int nodes2Triangles_host_vecs(int inode,
        HostSetInfoVecs& hostSetInfoVecs,
        CoordInfoVecs& coordInfoVecs,
        GeneralParams& generalParams,
        AuxVecs& auxVecs
        );
    

    int edge_swap_device_vecs (int iedge, 
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);

    int edge_swap_host_vecs (int iedge, 
        GeneralParams& generalParams,
        HostSetInfoVecs& hostSetInfoVecs,
        LinearSpringInfoVecs& linearSpringInfoVecs,
        BendingTriangleInfoVecs& bendingTriangleInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);
	
    void transferHtoD(GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);
    
    void transferDtoH(GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);

    void triangles2Triangles_host_vecs(
    int elem,
    HostSetInfoVecs& hostSetInfoVecs,
    CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams,
    AuxVecs& auxVecs);

    void gradient_weakening_update_host_vecs(double sigma,
        //double max_height_index,
        double max_height_x,
        double max_height_y,
        double max_height_z,
        double distance_to_boundary,
        double distance_to_boundary_max,
        GeneralParams& generalParams,
        CoordInfoVecs& coordInfoVecs,
        HostSetInfoVecs& hostSetInfoVecs);

    void gradient_weakening_update_host_vecs_tip(double sigma,
    //double max_height_index,
    double max_height_x,
    double max_height_y,
    double max_height_z,
    double distance_to_boundary,
    double distance_uniform_weak,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs);

    void chemSignal_weakening_update_host_vecs(double max_conc,
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    HostSetInfoVecs& hostSetInfoVecs);
    
    void make_ele_structure_PDE(ELEMENT element, int id, CoordInfoVecs& coordInfoVecs,
        AreaTriangleInfoVecs& areaTriangleInfoVecs);
    void make_rbc_structure(RBC* rbc, CoordInfoVecs& coordInfoVecs, GeneralParams& generalParams);
    //To use make_rbc_structure function, make sure you first initialize all the memory necessary to write the data structure.
    void adv_fw_LDG_update_u(std::list<ELEMENT>::iterator& it_ele,
                            std::list<ELEMENT>::iterator& n_it_ele,
                            std::list<ELEMENT>::iterator* it_adj_eles,
                            std::list<ELEMENT>::iterator* n_it_adj_eles,
                            RBC* rbc,
                            double dt,
                            double time,
                            double sqrt_diff_coef,
                            // double sink_coef,
                            // double source_coef,
                            double k_0,
                            double k_1,
                            double k_2,
                            double k_3,
                            double k_4,
                            double k_ss,
                            double beta,
                            double gamma,
                            double q1,
                            double h,
                            int rk_iter,
                            int INCLUDE_ADV);
    void adv_fw_LDG_update_q(std::list<ELEMENT>::iterator& it_ele,
                            std::list<ELEMENT>::iterator& n_it_ele,
                            std::list<ELEMENT>::iterator* it_adj_eles,
                            std::list<ELEMENT>::iterator* n_it_adj_eles,
                            RBC* rbc,
                            double dt,
                            double sqrt_diff_coeff,
                            int rk_iter);
    void Comp_ele_affine_trans(RBC* rbc, CoordInfoVecs& coordInfoVecs);
    void Comp_ele_mass_matrix(RBC* rbc, CoordInfoVecs& coordInfoVecs);
    void matrix_matrix_mult(double **mat,
                            double **matr,
                            int row,
                            int col,
                            double **ans);
    void matrix_vec_mult(double **mat,
                        double *vec,
                        int row,
                        int col,
                        double *ans);
    void soln_at_pt(double *dg_poly,
                    double *crds,
                    double *cent,
                    double sqrt_area,
                    double *soln);
    double inter_u_integr_7_quad_LDG(
                    double    *vt[], // vertices of the tri after affine mapping. 
                    double    *cent,  // after affine mapping
                    double    area, 
                    double    *soln_u, // soln polynomial: u of LDG
                    double    crds[][3],
                    double    fluxu[][4]);
    double inter_integr_7_quad_flux_LDG(
        double    *vt[], // vertices of the tri after affine mapping. 
        double    *cent,
        double    area,
        double    **affine_trans, 
        double    **inv_affine_tran,
        double    *q[],       //input soln: flux q of LDG in local coords.
        double    crds[][3],  //output. quadrature pts where vec{q} is evaluated. 
        double    fluxx[][4], //output: [][# eq]
        double    fluxy[][4],
	    double    fluxz[][4]);
    void inter_integr_7_quad_ver2(
	    double        *vt[], // vertices of the tri after affine mapping,
        double        *cent,
        double        area,
        int           indx,
        double        crds[][3],
        double        fluxx[][4], //[#quadrature][# eq]
        double        fluxy[][4],
        double        *ans);
    void inter_u_integr_7_quad_ver2(
	    double        *vt[], // vertices of the tri after affine mapping,
        double        *cent,
        double        area, 
        int           indx,
        double        crds[][3],
        double        u_flux[][4],
        double        *ans);
   
    void edge_integr_q_hat_flux(
        std::list<ELEMENT>::iterator &it_ele,
        std::list<ELEMENT>::iterator &it_adj_ele,
        double                  *prev_soln,        // prev rk-stepping u in element local crds
        double                  *prev_nb_soln,     // prev rk-stepping u in neighbor element local crds
        double                  *q_prev_soln[],    // prev rk-stepping q in element local crds
        double                  *q_prev_nb_soln[], // prev rk-stepping q in neighbor element local crds
        int                     side,
        int                     indx,
        double                  *e_q_int, // \int_{edge} q_hat\cdot n v_h ds
        // double                  u_hat[][10], // on edge u_hat, [# eq.][# quadra pt]
        double                  dt,
        int                     rk_iter,
        double                  *affine_pt[], // vertices of the tri after affine mapping
        double                  *aff_cent,
        double                  **adj_aff_time_inv_aff, 
                                 //neighbors' affine_tran times inverse of affine_tran of current ele (it_ele).
        double                  **aff_time_inv_adj_aff, 
                                 //affine_tran of current ele. times inverse of neighbors' affine_tran 
        double                  *adj_aff_cent, // center of neighbor element after affine trans. in
                                               // its own local coords.
        double                  **inv_affine_tran,
        double                  **inv_adj_affine_tran,
	    RBC                     *rbc);
    void edge_integr_new_u_hat_flux(
        std::list<ELEMENT>::iterator &it_ele,
        std::list<ELEMENT>::iterator &it_adj_ele,
        double                  *cur_soln, // in element local crds
        double                  *cur_nb_soln, // in neighbor element local crds
        int                     side,
        int                     indx,     
        double                  *e_u_int, // \int_{edge} u_hat n \cdot (v_h) ds
        double                  dt,
        int                     rk_iter,
        double                  *affine_pt[], // vertices of the tri after affine mapping
        double                  *aff_cent,
        double                  **adj_aff_time_inv_aff, 
                                 //neighbors' affine_tran times inverse of affine_tran of current
        double                  **aff_time_inv_adj_aff, 
                                 //affine_tran of current times inverse of neighbors' affine_tran 
        double                  *adj_aff_cent, // center of neighbor element after affine trans.
        double                  **inv_affine_tran,
        double                  **inv_adj_affine_tran,
	    RBC                     *rbc,
        bool                    *it_f_flag,
        bool                    *it_adj_f_flag,
        double                  it_hat_flux_store[][4][2],
        double                  it_adj_hat_flux_store[][4][2]);

    void inverse_matrix(
        double      **mat,
        int         size,
        double      **inv);

    void inverse_matrix_gj(
        double      **mat,
        int        size,
        double      **inv);
    void gaussj(
        long double **a,
        int         n,
        long double **b,
        int         m);
    void grad_vh_loc_basis(
        double  *crds,
        double  *cent,
        double  sqrt_area,
        int     indx,
        double  *val);
    void div_vec_vh_loc_basis(
        double  *crds,
        double  *cent,
        double  sqrt_area,
        int     indx,
        double  *val);
    void Tri_affine_trans(
	double   *pt0,
        double   *pt1,
        double   *pt2,
        double   *norm, 
        double   **affine_trans);
    double triangle_area_3d(
       double   *p1,
       double   *p2,
       double   *p3);
    void comp_mass_matrix(
        int      n_coeff,
        double   *pt0,
        double   *pt1,
        double   *pt2,
        double   *cent,
        double   sqrt_area,
        double   **mass_m);
    void comp_mass_matrix_of_vec_vh(
        int      n_coeff,
        double   *pt0,
        double   *pt1,
        double   *pt2,
        double   *cent,
        double   sqrt_area,
        double   **mass_m);
    double vh_val_loc_div_free_basis(
        double *crds,
        double *cent,
        double sqrt_area,
        int    indx);
    int LDG_C12_vector(
	double *n_pos,
	double *n_neg,
	int    id_pos,
	int    id_neg,
        int    dim, 
	double *C12);
    void vec_vh_val_loc_div_free_basis(
        double  *crds,
        double  *cent,
        double  sqrt_area,
        int     indx,
        double  *ans);
    double enorm0_3d(
       double   *p0,
       double   *p1);
       void tri_quadrature_16_pts(
        double       *pcrds0,
        double       *pcrds1,
        double       *pcrds2,
        double       crds[][2]);
    double normalized_val(
        double   crds[][2],
        double   *cent,
        double   sqrt_area,
        int      pos,
        int      indx);
    void  normalized_vec_vh_val(
        double  crds[][2],
        double  *cent,
        double  sqrt_area,
        int     pos,
        int     indx,
        double  *ans);
    void LDG_Surf_diffusion(
        RBC *rbc, //will throw in the address of the rbc struct
        RBC *n_rbc,
        double sqrt_D,
        double k_0,
        double k_1,
        double k_2,
        double k_3,
        double k_4,
        double k_ss,
        double beta,
        double gamma,
        double q1,
        double h,
        // double sink_coef,
        // double source_coef,
        double fr_time,
        double time_step_factor,
        //int        *iperm,
        //double     *dh,
        //Wave       *wv,
        //Wave       *nwv,
        //Front      *fr,
        //Front      *nfr,
        double     dt,
        int        INCL_ADV
        );
    //void Allocate_memory(RBC *rbc);
    double current_area(RBC *rbc);
    void build_ele_adj_by_iter(RBC* rbc, CoordInfoVecs& coordInfoVecs);
    void comput_aff_conormals(RBC* rbc);
    void comput_aff_cent_vtx(RBC* rbc);
    
    // double derivative_P4_surf(
	// double   *surf,
	// double   *x,
    //     int      indx);
   
    void accurate_sphere_diffusion_initializer(
        RBC                      *rbc,
        std::list<ELEMENT>::iterator  &it_ele, 
        double                   *dg_rho_soln,
        double                   *dg_q[]);
    
    // void Init_cell_intfc_vel(
	// RBC* rbc
    // //Wave     *wv,
    //   //  Front    *fr
    //   );
    double accurate_diffusion_sphere(
	double                   *glb_cell_cent, 
	double                   *glb_crds,
        double                   avg_radi,
        std::list<ELEMENT>::iterator  &it_ele,
	double                   time);
    void LDG_Error_of_surf_diffusion_test(
	RBC      *rbc, 
	double   time, 
	double   *diff_error);
    void LDG_Error_of_surf_diffusion_test_on_ele(
	RBC      *rbc, 
	std::list<ELEMENT>::iterator it_ele,
	double   time);

    void inter_integr_7_quad_ver3(
    double        *prev_soln,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans);

    void inter_integr_7_quad_ver4(
    double        source_coef,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans);

    void inter_integr_7_quad_ver5(
    double        k_0,
    double        beta,
    double        q1,
    double        u,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans);

    void inter_integr_7_quad_ver6(
    double        k_1,
    double        beta,
    double        gamma,
    double        q1,
    double        h,
    double        u,
    double        *prev_soln,
    double        *cent,
    double        area,
    int           indx,
	double        crds[][3],
	double        *ans);

    void LDG_Surface_Diffusion_Solve(
    double edgeswap_iteration,
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc);
    void LDG_Surface_Diffusion_Initialize(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc,
    uint mem_prealloc
    );
    void LDG_Surface_Diffusion_Structure_Rebuild(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc
    );
    void LDG_Surface_Diffusion_Structure_Rebuild_postGrowth(
    CoordInfoVecs& coordInfoVecs,
    GeneralParams& generalParams,
    HostSetInfoVecs& hostSetInfoVecs,
    AuxVecs& auxVecs,
	RBC      *rbc,
    RBC      *n_rbc
    );
};



#endif
