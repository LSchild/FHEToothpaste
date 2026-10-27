//
// Created by Leonard on 7/23/26.
//

#include "backend/backend.h"
#include "backend/hexl_worker.h"
#include "backend/trivial_worker.h"
#include "backend/solinas_worker.h"

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

    auto worker = std::make_shared<HexlWorker>(modulus, dimension);
    WORKER_MAP.emplace(stupid_hash, worker);

    return worker;
}

std::shared_ptr<MathWorker> SelectWorkerArbitraryDimension(uint64_t modulus) {

    const uint64_t min_dim_bits = 4;
    const uint64_t max_dim_bits = 17;

    for(uint64_t i = min_dim_bits; i < max_dim_bits; i++) {

        auto tmp_dim = 1ull << i;
        uint64_t stupid_hash = modulus ^ (tmp_dim | (~tmp_dim << 32));
        if (WORKER_MAP.contains(stupid_hash)) {
            return WORKER_MAP[stupid_hash];
        }
    }

    const uint64_t dim = 1024;
    uint64_t stupid_hash = modulus ^ (dim | (dim << 32));

    if (modulus == (-(1ull << 32) + 1)) {
        auto solinas_worker = std::make_shared<SolinasBE::SolinasWorker>(dim);
        WORKER_MAP[stupid_hash] = solinas_worker;
        return solinas_worker;
    }

    auto worker = std::make_shared<TrivialWorker>(modulus, 0);
    WORKER_MAP.emplace(stupid_hash, worker);

    return worker;

}