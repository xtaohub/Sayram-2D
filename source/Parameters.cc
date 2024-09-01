#include <iostream>
#include <filesystem>
#include "Parameters.h"
#include "Ini_reader.h"

namespace fs = std::filesystem; 

Parameters::Parameters(int argc, char** argv){

  handle_main_input(argc, argv);
  read_inp_file(); 

  output_path_ = "./output/" + run_id() + "/"; 
  fs::create_directories(output_path_); 

 // copy the parameter file 
  string paras_file = "./output/" + run_id() + "/" + run_id() + ".ini";
  string command = "cp " + inp_file() + " " +  paras_file;
  system(command.c_str());
}

void Parameters::handle_main_input(int argc, char* argv[]){
  switch (argc) {
  case 1:
    inp_file_ = "p.ini"; 
    break;
    
  case 2:
    inp_file_.assign(argv[1]); 
    break;
    
  default:
    std::cerr << "NParas() > 2! This program takes at most one argument: the parameter file name." << std::endl; 
    exit(1); 
  }
}

void Parameters::read_inp_file(){

  Ini_reader ireader(inp_file()); 

  ireader.set_section("basic");

  ireader.read("run_id", &run_id_);
  ireader.read("nalpha0", &nalpha0_);
  ireader.read("nE", &nE_);

  assert(nE_ > 0 && nalpha0_ > 0);
   
  ireader.read("L", &L_);
  ireader.read("alpha0_min_bct", &alpha0_min_bct_); 

  alpha0_lc_ = asin(pow(pow(L_,5)*(4*L_-3), -0.25)); 

  if (alpha0_min_bct_ == 0) 
    alpha0_min_ = 0.0; 
  else 
    alpha0_min_ = alpha0_lc_; 

  alpha0_max_ = gPI/2.0; 

  ireader.read("Emin", &Emin_);
  ireader.read("Emax", &Emax_);

  pmin_ = e2p(Emin_, gE0);
  pmax_ = e2p(Emax_, gE0);

  ireader.read("T", &T_);
  ireader.read("nsteps", &nsteps_);

  ireader.set_section("diagnostics");

  ireader.read("nplots", &nplots_); 
  save_every_step_ = nsteps_ / nplots_; 
  nsteps_ = save_every_step_ * nplots_; 

  ireader.set_section("diffusion_coefficients"); 

  ireader.read("dID", &dID_);
  ireader.read("nalpha0_D", &nalpha0_D_);
  ireader.read("alpha0_min_D", &alpha0_min_D_);
  ireader.read("alpha0_max_D", &alpha0_max_D_);

  alpha0_min_D_ = alpha0_min_D_ * gPI / 180.0; 
  alpha0_max_D_ = alpha0_max_D_ * gPI / 180.0; 
  dalpha0_D_ = (alpha0_max_D_ - alpha0_min_D_)/(nalpha0_D_ - 1); 

  ireader.read("nE_D", &nE_D_);
  ireader.read("Emin_D", &Emin_D_);
  ireader.read("Emax_D", &Emax_D_);

  dlogE_D_ = (log(Emax_D_) - log(Emin_D_)) / (nE_D_ - 1); 

}
