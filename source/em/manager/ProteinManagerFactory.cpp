#include <em/manager/ProteinManagerFactory.h>
#include <em/manager/SmartProteinManager.h>
#include <em/manager/SimpleProteinManager.h>
#include <data/Molecule.h>
#include <settings/EMSettings.h>

std::unique_ptr<em::managers::ProteinManager> em::factory::create_manager(const ImageStackBase* images) {
    switch (settings::em::fixed_weights) {
        case true: return std::make_unique<em::managers::SimpleProteinManager>(images);
        case false: return std::make_unique<em::managers::SmartProteinManager>(images);
    }
}