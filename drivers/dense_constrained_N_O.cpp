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

struct PairSort
{
    GQCP::Eigenpair pair;
    size_t origin;
    double mul;
    double entropy;

    PairSort(GQCP::Eigenpair pair, size_t origin, double mul, double entropy) : pair(std::move(pair)), origin(origin), mul(mul), entropy(entropy) {};
    PairSort(const PairSort& pairs) :  pair(pairs.pair), origin(pairs.origin), mul(pairs.mul), entropy(pairs.entropy) {};


    bool operator<(const PairSort& x) const {
        if (this->pair.isEqual(x.pair)){
            return false;
        }
        return this->pair.get_eigenvalue() < x.pair.get_eigenvalue();
    };

    bool operator>(const PairSort& x) const {
        if (this->pair.isEqual(x.pair)){
            return false;
        }
        return this->pair.get_eigenvalue() > x.pair.get_eigenvalue();
    };

    bool operator==(const PairSort& x) const {
        return (this->pair.isEqual(x.pair));
    };

};


int main(int argc, char** argv) {

    /**
     *  LOGISTICS
     */
    // Time the EXE
    auto start = std::chrono::high_resolution_clock::now();
    size_t target = 5;
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

    std::ostringstream distance_string_precursor;
    distance_string_precursor << std::setprecision(2) << distance;
    std::string distance_string = distance_string_precursor.str();

    std::vector<std::ofstream> outputfiles;

    for (size_t i = 0; i < target; i++) {
        std::string output_filename = atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_dense_" + basisset + name + ".out" + std::to_string(i);
        std::ofstream output_file;
        output_file.open(output_filename, std::fstream::out);
        outputfiles.push_back(std::move(output_file));
    }
    std::string output_filename_log = atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_dense_" + basisset + name + ".log" ;
    std::ofstream output_log;
    output_log.open(output_filename_log, std::fstream::out);
    output_log << "init" <<std::endl;
    output_log << "Version: " << std::setprecision(15) << GQCP_GIT_SHA1 <<  std::endl;
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

    GQCP::AtomicDecompositionParameters adp(molecule, basisset);

    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(molecule, basisset);  // in the AO basis
    auto K = mol_ham_par.get_K();
    mol_ham_par.LowdinOrthonormalize();

    // chose first half of AO as constrain targets
    std::vector<size_t> AOlist;
    for (size_t i = 0; i < K/2; i++) {
        AOlist.push_back(i);
    }

    Eigen::Map<Eigen::Matrix<size_t, Eigen::Dynamic, 1>> bfmap (AOlist.data(), AOlist.size());
    GQCP::VectorXs bfsv (bfmap);

    output_log << "selected BF: " << std::setprecision(15) << std::endl << bfsv.transpose() << std::endl;

    GQCP::FrozenProductFockSpace fock_space (K, N_alpha, N_beta, X);
    GQCP::DenseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = fock_space.get_dimension(); // request to store all eigenpairs

    /**
     *  CALCULATIONS
     */
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
    GQCP::FrozenCoreFCI fci (fock_space);
    GQCP::RDMCalculator rdm_calculator (fock_space);
    auto mulliken_operator = mol_ham_par.calculateMullikenOperator(AOlist);

    for (size_t i = 0; i < lambdas.rows(); i++) { // lambda iterations

        std::string lambdastr = "_lambda" + std::to_string(lambdas(i));
        std::string dense_filename;
        std::ofstream dense_log;
        dense_filename =
                atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_dense_" + basisset + lambdastr + name + ".full"; // NOLINT(performance-inefficient-string-concatenation)
        dense_log.open(dense_filename, std::fstream::out);

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
        std::cout << "TOTAL DENSE SOLVE TIME" << " : " << seconds1 << " seconds" << std::endl;

        const auto& all_pairs = ci_solver.get_eigenpairs();
        std::vector<PairSort> sorting_vector;
        sorting_vector.reserve(fock_space.get_dimension());

        double total = 0;
        for (size_t j = 0; j < fock_space.get_dimension(); j++) {

            const GQCP::Eigenpair& pair = all_pairs[j];

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
            total += seconds2;

            double mul = calculateExpectationValue(mulliken_operator, D);
            GQCP::WaveFunction wavefunction (fock_space, fci_coefficients);
            double entropy = wavefunction.calculateShannonEntropy();
            double fci_energy = en + internuclear_repulsion_energy + lambdas(i) * mul;

            const auto& T = mol_ham_par.get_T_total();
            auto D_ao = T * D * T.adjoint();
            d.fourModeMultiplication<double>(T.adjoint().conjugate(), T.adjoint(), T.adjoint().conjugate(), T.adjoint());

            GQCP::Eigenpair constrained_pair(fci_energy, fci_coefficients);
            sorting_vector.emplace_back(constrained_pair, j, mul, entropy);
        }
        std::cout << "RDM TIME" << " : " << total << " seconds" << std::endl;
        auto start3 = std::chrono::high_resolution_clock::now();

        std::sort(sorting_vector.begin(), sorting_vector.end(), []( const PairSort &left, const PairSort &right )
        { return ( left.pair.get_eigenvalue() < right.pair.get_eigenvalue() ); } );

        auto stop3 = std::chrono::high_resolution_clock::now();
        // Process the chrono time and output
        auto elapsed_time3 = stop3 - start3;           // in nanoseconds
        auto seconds3 = elapsed_time3.count() / 1e9;  // in seconds
        std::cout << "sort time" << " : " << seconds3 << "seconds" << std::endl;

        size_t counter = 0;
        for (auto const& sorted : sorting_vector) {

            if (counter < target) {

                const auto& fci_coefficients = sorted.pair.get_eigenvector();

                rdm_calculator.set_coefficients(fci_coefficients);
                GQCP::OneRDM<double> D = rdm_calculator.calculate1RDMs().one_rdm;
                GQCP::TwoRDM<double> d = rdm_calculator.calculate2RDMs().two_rdm;

                const auto& T = mol_ham_par.get_T_total();
                auto D_ao = T * D * T.adjoint();
                d.fourModeMultiplication<double>(T.adjoint().conjugate(), T.adjoint(), T.adjoint().conjugate(), T.adjoint());

                double en_A = GQCP::calculateExpectationValue(adp.fragment_parameters[0], D_ao, d);
                double en_B = GQCP::calculateExpectationValue(adp.fragment_parameters[1], D_ao, d);
                double en_AA = GQCP::calculateExpectationValue(adp.net_atomic_parameters[0], D_ao, d);
                double en_BB = GQCP::calculateExpectationValue(adp.net_atomic_parameters[1], D_ao, d);
                double en_AB = GQCP::calculateExpectationValue(adp.interaction_parameters[0], D_ao, d);



                auto& output_file = outputfiles[counter];
                output_file << std::setprecision(15) << sorted.pair.get_eigenvalue()  << "\t" << lambdas(i) << "\t" << sorted.mul << "\t" << sorted.entropy << "\t"
                    << en_A << "\t" << en_AA << "\t" << en_B << "\t" << en_BB << "\t" << en_AB << "\t"
                        << sorted.origin << "\t" << counter << std::endl;
            }

            dense_log << std::setprecision(15) << sorted.pair.get_eigenvalue()  << "\t" << lambdas(i) << "\t" << sorted.mul << "\t" << sorted.entropy << "\t"
                    << "\t"
                    << sorted.origin << "\t" << counter << std::endl;
            counter++;

            if (counter == 100) {
                break;
            }
        }

        dense_log.close();
    }

    output_log << "mullikenoperator: " << std::setprecision(15) << std::endl << mulliken_operator << std::endl;
    output_log << "Total C: " << std::setprecision(15) << std::endl << mol_ham_par.get_T_total() << std::endl;
    output_log << "Basis set used: " << std::setprecision(15) << basisset << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    // Process the chrono time and output
    auto elapsed_time = stop - start;           // in nanoseconds
    auto seconds = elapsed_time.count() / 1e9;  // in seconds
    output_log << "TOTAL EXECUTABLE TIME" << " : " << seconds << " seconds" << std::endl;

    for (auto& output_file : outputfiles) {
        output_file.close();
    }
    output_log.close();
    return 0;
}