//
// Created by Leonard on 7/23/26.
//

#include "backend/backend.h"
#include "backend/hexl_worker.h"
#include "backend/solinas.h"

#include <unordered_map>

static std::unordered_map<uint64_t, std::shared_ptr<MathWorker>> WORKER_MAP = {

};

std::shared_ptr<MathWorker> SelectWorker(uint64_t modulus, uint64_t dimension) {

    // TODO: use better hashing function
    uint64_t stupid_hash = modulus ^ (dimension | (~dimension << 32));
    if (WORKER_MAP.contains(stupid_hash)) {
        return WORKER_MAP[stupid_hash];
    }

    if (modulus == (-(1ull << 32) + 1)) {
        auto solinas_worker = std::make_shared<SolinasBE::SolinasWorker>(dimension);
        WORKER_MAP[stupid_hash] = solinas_worker;
        return solinas_worker;
    }

    // TODO Future: select worker based on Q & dimensio value
    auto worker = std::make_shared<HexlWorker>(modulus, dimension);
    WORKER_MAP.emplace(stupid_hash, worker);

    return worker;
}