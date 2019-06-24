#include <utility>

/**
 *  This executable calculates NO+ for a given basis set, at a given intra-molecular distance and constraining a given range of lambdas
 *  Eigenvalue problem is solved Dense, all eigenvalues are logged.
 */


#include <fstream>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <gqcp.hpp>
#include <HamiltonianParameters/AtomicDecompositionParameters.hpp>
#include <HamiltonianBuilder/FrozenCoreFCI.hpp>
#include <algorithm>

namespace po = boost::program_options;


int main(int argc, char** argv) {

    /**
     *  LOGISTICS
     */
    // Time the EXE
    auto start = std::chrono::high_resolution_clock::now();
    // Molecule specifications
    std::string atom_str1 = "N";
    std::string atom_str2 = "O";
    size_t N_alpha = 7;
    size_t N_beta = 7;

    // Frozencores
    size_t X = 0;

    // Output parsing
    std::string name;
    std::time_t now = std::time(0);
    std::string total_tag = std::to_string(now);

    // Input processing
    std::string basisset;
    double distance;
    std::string lambda_string;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
                ("help,h", "print help messages")
                ("distance,d", po::value<double>(&distance)->required(), "intranuclear distance")
                ("constraint,c", po::value<std::string>(&lambda_string)->required(), "cs of all lambdas")
                ("basis,s", po::value<std::string>(&basisset)->required(), "name of the basis set")
                ("frozencores,x", po::value<size_t>(&X)->default_value(0), "freeze amount of orbitals")
                ("name,e", po::value<std::string>(&name)->default_value(""), "name extension for the file");

        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }

    std::vector<std::string> splitted_line_lambda;
    boost::split(splitted_line_lambda, lambda_string, boost::is_any_of(","));

    Eigen::VectorXd lambdas(splitted_line_lambda.size());
    size_t index = 0;
    for (const std::string& x : splitted_line_lambda) {
        lambdas(index) = std::stod(x);
        index++;
    }

    // Create and open a file

    std::string output_filename_log = name + ".log";
    std::string output_filename_bin = name + ".out";
    std::ofstream output_log;
    std::ofstream output_file;
    output_file.open(output_filename_bin, std::fstream::out);
    output_log.open(output_filename_log, std::fstream::out);
    output_log << "init" <<std::endl;
    output_log << "Version: " << std::setprecision(15) << "placeholdergit" <<  std::endl;
    output_log << "Frozencore? : " << std::setprecision(15) << std::endl << X << std::endl;
    output_log << "selected lambdas: " << std::setprecision(15) << std::endl << lambdas.transpose() << std::endl;

    /**
     *  ENVIRONMENT
     */

    std::vector<GQCP::Atom> atom_list;

    GQCP::Atom atom1 (GQCP::elements::elementToAtomicNumber(atom_str1), -distance/2, 0, 0);
    GQCP::Atom atom2 (GQCP::elements::elementToAtomicNumber(atom_str2), distance/2, 0, 0);

    atom_list.push_back(atom1);
    atom_list.push_back(atom2);

    // +1 charge we are hard coding NO+
    GQCP::Molecule molecule (atom_list, +1);

    GQCP::AtomicDecompositionParameters adp = GQCP::AtomicDecompositionParameters::Nuclear(molecule, basisset);

    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(molecule, basisset);  // in the AO basis
    auto K = mol_ham_par.get_K();

    try {
        // define individual atoms as molecular fractions
        GQCP::Molecule mol_fraction1(std::vector<GQCP::Atom>{atom1}, +1);
        GQCP::Molecule mol_fraction2(std::vector<GQCP::Atom>{atom2});

        auto ham_par1 = GQCP::HamiltonianParameters<double>::Molecular(mol_fraction1, basisset);
        auto ham_par2 = GQCP::HamiltonianParameters<double>::Molecular(mol_fraction2, basisset);

        // Perform DIIS RHF for individual fractions
        GQCP::DIISRHFSCFSolver diis_scf_solver1 (ham_par1, mol_fraction1, 6, 6, 1e-12, 500);
        GQCP::DIISRHFSCFSolver diis_scf_solver2 (ham_par2, mol_fraction2, 6, 6, 1e-12, 500);
        diis_scf_solver1.solve();
        diis_scf_solver2.solve();
        auto rhf1 = diis_scf_solver1.get_solution();
        auto rhf2 = diis_scf_solver2.get_solution();

        // Retrieve transformation from the solutions and transform the Hamiltonian parameters
        size_t K1 = ham_par1.get_K();
        size_t K2 = ham_par2.get_K();

        GQCP::SquareMatrix<double> T = Eigen::MatrixXd::Zero(K, K);
        T.topLeftCorner(K1, K1) += rhf1.get_C();
        T.bottomRightCorner(K2, K2) += rhf2.get_C();
        mol_ham_par.basisTransform(T);

        // Perform DIIS with the new basis if this fails Lodwin orthonormalize
        try {
            GQCP::DIISRHFSCFSolver diis_scf_solver (mol_ham_par, molecule, 6, 6, 1e-12, 500);
            diis_scf_solver.solve();
            auto rhf = diis_scf_solver.get_solution();
            mol_ham_par.basisTransform(rhf.get_C());

        } catch (const std::exception& e) {
            output_log << "Lodwin Orthonormalized" << std::endl;
            mol_ham_par.LowdinOrthonormalize();
        }

    } catch (const std::exception& e) {
        output_log << e.what() << std::endl;
        output_file.close();
        output_log.close();
        return 1;
    }

    // chose first half of AO as constrain targets
    std::vector<size_t> AOlist;
    for (size_t i = 0; i < K/2; i++) {
        AOlist.push_back(i);
    }

    Eigen::Map<Eigen::Matrix<size_t, Eigen::Dynamic, 1>> bfmap (AOlist.data(), AOlist.size());
    GQCP::VectorXs bfsv (bfmap);

    output_log << "selected BF: " << std::setprecision(15) << std::endl << bfsv.transpose() << std::endl;

    GQCP::FrozenProductFockSpace fock_space (K, N_alpha, N_beta, X);
    GQCP::DavidsonSolverOptions solver_options (fock_space.HartreeFockExpansion());

    /**
     *  CALCULATIONS
     */
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
    GQCP::FrozenCoreFCI fci (fock_space);
    GQCP::RDMCalculator rdm_calculator (fock_space);
    auto mulliken_operator = mol_ham_par.calculateMullikenOperator(AOlist);

    for (size_t i = 0; i < lambdas.rows(); i++) { // lambda iterations
        // Created constrained ham_par
        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator, lambdas(i));

        GQCP::CISolver ci_solver (fci, constrained_ham_par);

        auto start1 = std::chrono::high_resolution_clock::now();
        // SOLVE
        try {
            ci_solver.solve(solver_options);
        } catch (const std::exception& e) {
            output_log << e.what() << "lambda: " << lambdas(i);
            continue;
        }

        auto stop1 = std::chrono::high_resolution_clock::now();
        // Process the chrono time and output
        auto elapsed_time1 = stop1 - start1;           // in nanoseconds
        auto seconds1 = elapsed_time1.count() / 1e9;  // in seconds
        std::cout << "TOTAL DAVIDSON SOLVE TIME" << " : " << seconds1 << " seconds" << std::endl;

        const GQCP::Eigenpair& pair = ci_solver.get_eigenpair();


        auto start2 = std::chrono::high_resolution_clock::now();

        const auto& fci_coefficients = pair.get_eigenvector();
        double en = pair.get_eigenvalue();
        rdm_calculator.set_coefficients(fci_coefficients);
        GQCP::OneRDM<double> D = rdm_calculator.calculate1RDMs().one_rdm;
        GQCP::TwoRDM<double> d = rdm_calculator.calculate2RDMs().two_rdm;

        auto stop2 = std::chrono::high_resolution_clock::now();
        // Process the chrono time and output
        auto elapsed_time2 = stop2 - start2;           // in nanoseconds
        auto seconds2 = elapsed_time2.count() / 1e9;  // in seconds
        std::cout << "RDM TIME" << " : " << seconds2 << " seconds" << std::endl;

        double mul = calculateExpectationValue(mulliken_operator, D);
        GQCP::WaveFunction wavefunction (fock_space, fci_coefficients);
        double entropy = wavefunction.calculateShannonEntropy();
        double fci_energy = en + internuclear_repulsion_energy + lambdas(i) * mul;

        const auto& T = mol_ham_par.get_T_total();
        D.basisTransform<double>(T.adjoint());;
        d.basisTransform<double>(T.adjoint());

        double en_A = GQCP::calculateExpectationValue(adp.get_atomic_parameters()[0], D, d);
        double en_B = GQCP::calculateExpectationValue(adp.get_atomic_parameters()[1], D, d);
        double en_AA = GQCP::calculateExpectationValue(adp.get_net_atomic_parameters()[0], D, d);
        double en_BB = GQCP::calculateExpectationValue(adp.get_net_atomic_parameters()[1], D, d);
        double en_AB = GQCP::calculateExpectationValue(adp.get_interaction_parameters()[0], D, d);


        output_file << std::setprecision(15) << fci_energy << "\t" << lambdas(i) << "\t" << mul << "\t" << entropy << "\t"
                    << en_A << "\t" << en_AA << "\t" << en_B << "\t" << en_BB << "\t" << en_AB
                    << std::endl;

        solver_options.X_0 = fci_coefficients;
    }

    output_log << "mullikenoperator: " << std::setprecision(15) << std::endl << mulliken_operator << std::endl;
    output_log << "Total C: " << std::setprecision(15) << std::endl << mol_ham_par.get_T_total() << std::endl;
    output_log << "Basis set used: " << std::setprecision(15) << basisset << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    // Process the chrono time and output
    auto elapsed_time = stop - start;           // in nanoseconds
    auto seconds = elapsed_time.count() / 1e9;  // in seconds
    output_log << "TOTAL EXECUTABLE TIME" << " : " << seconds << " seconds" << std::endl;

    output_file.close();
    output_log.close();

    return 0;
}