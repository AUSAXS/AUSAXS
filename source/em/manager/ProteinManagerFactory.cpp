#include <em/manager/ProteinManagerFactory.h>
#include <em/manager/SmartProteinManager.h>
#include <em/manager/SimpleProteinManager.h>
#include <data/Molecule.h>
#include <settings/EMSettings.h>

std::unique_ptr<em::managers::ProteinManager> em::factory::create_manager(observer_ptr<const ImageStackBase> images) {
    if (settings::em::fixed_weights) {
        return std::make_unique<em::managers::SimpleProteinManager>(images);
    } else {
        return std::make_unique<em::managers::SmartProteinManager>(images);
    }
}