#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <cmath>
#include "utils.h"
#include "common.h"

class Parameters{
public:
  Parameters(int argc, char** argv); 

  // ------- READ FROM INP FILE -------
  const string& run_id() const { return run_id_; }
  const string& inp_file() const { return inp_file_; }

  int nalpha0() const { return nalpha0_; }
  int nE() const { return nE_; }
  double alpha0_lc() const { return alpha0_lc_ * gPI / 180.0; }
  double alpha0_max() const { return alpha0_max_ * gPI / 180.0; }

  double Emin() const { return Emin_; }
  double Emax() const { return Emax_; }

  double T() const { return T_; }
  int nsteps() const { return nsteps_; }
  double dt() const { return T_ / nsteps_;}

  double pmin() const { return pmin_; }
  double pmax() const { return pmax_; }

  int nplots() const { return nplots_; }
  int save_every_step() const { return save_every_step_; }
  const string& output_path() const { return output_path_; }

  const string& dID() const { return dID_; }
  int nalpha0_D() const { return nalpha0_D_; }
  double alpha0_min_D() const { return alpha0_min_D_; }
  double alpha0_max_D() const { return alpha0_max_D_; }
  double dalpha0_D() const { return dalpha0_D_; }

  int nE_D() const { return nE_D_; }
  double Emin_D() const { return Emin_D_; }
  double Emax_D() const { return Emax_D_; }
  double dlogE_D() const { return dlogE_D_; }

private:
  string inp_file_; 

  string run_id_;
  int nalpha0_;
  int nE_;

  double alpha0_lc_;
  double alpha0_max_; 

  double Emin_;
  double Emax_; 

  double pmin_;
  double pmax_; 

  double T_;
  double nsteps_;
  double dt_; 

  int nplots_;
  int save_every_step_; 
  string output_path_; 

  string dID_;

  int nalpha0_D_;
  double alpha0_min_D_;
  double alpha0_max_D_;
  double dalpha0_D_; 

  int nE_D_;
  double Emin_D_;
  double Emax_D_;
  double dlogE_D_;  

  void handle_main_input(int argc, char* argv[]);
  void read_inp_file(); 
};

#endif /* PARAMETERS_H_ */
