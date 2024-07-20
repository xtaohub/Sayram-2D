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
  ireader.read("nalpha", &nalpha_);
  ireader.read("nE", &nE_);
  ireader.read("alpha_lc", &alpha_lc_);
  ireader.read("alpha_max", &alpha_max_);
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
  ireader.read("nalpha_D", &nalpha_D_);
  ireader.read("alpha_min_D", &alpha_min_D_);
  ireader.read("alpha_max_D", &alpha_max_D_);
  ireader.read("nE_D", &nE_D_);
  ireader.read("Emin_D", &Emin_D_);
  ireader.read("Emax_D", &Emax_D_);

}
