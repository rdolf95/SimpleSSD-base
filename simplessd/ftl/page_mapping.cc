/*
 * Copyright (C) 2017 CAMELab
 *
 * This file is part of SimpleSSD.
 *
 * SimpleSSD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SimpleSSD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SimpleSSD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ftl/page_mapping.hh"

#include <algorithm>
#include <limits>
#include <random>

#include "util/algorithm.hh"
#include "util/bitset.hh"

SimpleSSD::Event refreshEvent;

namespace SimpleSSD {

namespace FTL {

PageMapping::PageMapping(ConfigReader &c, Parameter &p, PAL::PAL *l,
                         DRAM::AbstractDRAM *d)
    : AbstractFTL(p, l, d),
      pPAL(l),
      conf(c),
      lastFreeBlock(param.pageCountToMaxPerf),
      lastFreeBlockIOMap(param.ioUnitInPage),
      bReclaimMore(false),
      lastRefreshed(0) {
  blocks.reserve(param.totalPhysicalBlocks);
  table.reserve(param.totalLogicalBlocks * param.pagesInBlock);

  uint32_t initEraseCount = conf.readUint(CONFIG_FTL, FTL_INITIAL_ERASE_COUNT);

  for (uint32_t i = 0; i < param.totalPhysicalBlocks; i++) {
    //freeBlocks.emplace_back(Block(i, param.pagesInBlock, param.ioUnitInPage));
    freeBlocks.emplace_back(Block(i, param.pagesInBlock, param.ioUnitInPage, initEraseCount));
  }

  nFreeBlocks = param.totalPhysicalBlocks;

  status.totalLogicalPages = param.totalLogicalBlocks * param.pagesInBlock;

  // Allocate free blocks
  for (uint32_t i = 0; i < param.pageCountToMaxPerf; i++) {
    lastFreeBlock.at(i) = getFreeBlock(i);
  }

  lastFreeBlockIndex = 0;

  memset(&stat, 0, sizeof(stat));

  bRandomTweak = conf.readBoolean(CONFIG_FTL, FTL_USE_RANDOM_IO_TWEAK);
  bitsetSize = bRandomTweak ? param.ioUnitInPage : 1;

  float tmp = conf.readFloat(CONFIG_FTL, FTL_TEMPERATURE);
  float Ea = 1.1;
  float epsilon = conf.readFloat(CONFIG_FTL, FTL_EPSILON);
  float alpha = conf.readFloat(CONFIG_FTL, FTL_ALPHA);
  float beta = conf.readFloat(CONFIG_FTL, FTL_BETA);
  float gamma = conf.readFloat(CONFIG_FTL, FTL_GAMMA);
  float kTerm = conf.readFloat(CONFIG_FTL, FTL_KTERM);
  float mTerm = conf.readFloat(CONFIG_FTL, FTL_MTERM);
  float nTerm = conf.readFloat(CONFIG_FTL, FTL_NTERM);
  float sigma = conf.readFloat(CONFIG_FTL, FTL_ERROR_SIGMA);
  uint32_t seed = conf.readUint(CONFIG_FTL, FTL_RANDOM_SEED);


  errorModel = ErrorModeling(tmp, Ea, epsilon, alpha, beta, gamma,
                             kTerm, mTerm, nTerm, 
                             sigma, param.pageSize, seed);
}

PageMapping::~PageMapping() {
  //refreshStatFile.close();
}

bool PageMapping::initialize() {
  uint64_t nPagesToWarmup;
  uint64_t nPagesToInvalidate;
  uint64_t nTotalLogicalPages;
  uint64_t maxPagesBeforeGC;
  uint64_t tick;
  uint64_t valid;
  uint64_t invalid;
  FILLING_MODE mode;

  Request req(param.ioUnitInPage);

  debugprint(LOG_FTL_PAGE_MAPPING, "Initialization started");

  nTotalLogicalPages = param.totalLogicalBlocks * param.pagesInBlock;
  nPagesToWarmup =
      nTotalLogicalPages * conf.readFloat(CONFIG_FTL, FTL_FILL_RATIO);
  nPagesToInvalidate =
      nTotalLogicalPages * conf.readFloat(CONFIG_FTL, FTL_INVALID_PAGE_RATIO);
  mode = (FILLING_MODE)conf.readUint(CONFIG_FTL, FTL_FILLING_MODE);
  maxPagesBeforeGC =
      param.pagesInBlock *
      (param.totalPhysicalBlocks *
           (1 - conf.readFloat(CONFIG_FTL, FTL_GC_THRESHOLD_RATIO)) -
       param.pageCountToMaxPerf);  // # free blocks to maintain

  if (nPagesToWarmup + nPagesToInvalidate > maxPagesBeforeGC) {
    warn("ftl: Too high filling ratio. Adjusting invalidPageRatio.");
    nPagesToInvalidate = maxPagesBeforeGC - nPagesToWarmup;
  }

  debugprint(LOG_FTL_PAGE_MAPPING, "Total logical pages: %" PRIu64,
             nTotalLogicalPages);
  debugprint(LOG_FTL_PAGE_MAPPING,
             "Total logical pages to fill: %" PRIu64 " (%.2f %%)",
             nPagesToWarmup, nPagesToWarmup * 100.f / nTotalLogicalPages);
  debugprint(LOG_FTL_PAGE_MAPPING,
             "Total invalidated pages to create: %" PRIu64 " (%.2f %%)",
             nPagesToInvalidate,
             nPagesToInvalidate * 100.f / nTotalLogicalPages);

  req.ioFlag.set();

    //setup refresh
  //uint64_t random_seed = conf.readUint(CONFIG_FTL, FTL_RANDOM_SEED);
  uint32_t num_queue = 200;
  //uint32_t filter_size = conf.readUint(CONFIG_FTL, FTL_REFRESH_FILTER_SIZE);
  //debugprint(LOG_FTL_PAGE_MAPPING, "Refresh setting start. The number of bloom filters: %u", num_bf);
  //debugprint(LOG_FTL_PAGE_MAPPING, "Refresh threshold error count: %u", param.pageSize / 1000);
  
  //deque version
  for (uint32_t i=0; i<num_queue; i++){
    refreshQueues.push_back(std::deque<uint32_t>());
    checkedQueues.push_back(std::deque<uint32_t>());
  }
  insertedBlockCheck = Bitset(param.totalPhysicalBlocks); 
  // TODO : 8 = # of pages for each layer
  insertedBlockCheck.reset();

  // total physical blocks : total blocks in SSD
  // total logical blocks : total blocks that host can use (smaller than physical because of OP)

  debugprint(LOG_FTL_PAGE_MAPPING, "DIVCEIL(param.totalLogicalBlocks * param.pagesInBlock, 8): %u", DIVCEIL(param.totalLogicalBlocks * param.pagesInBlock, 8));
  debugprint(LOG_FTL_PAGE_MAPPING, "insertedLayerCheck.size(): %u", insertedBlockCheck.size());

  

  // set up periodic refresh event
  refresh_period = 3600000000000000*1;  // 1hour *24 = 1day
  
  if (refresh_period > 0) {
    refreshEvent = engine.allocateEvent([this](uint64_t tick) {
      refresh_event(tick);

      /*engine.scheduleEvent(
          refreshEvent,
          tick + conf.readUint(CONFIG_FTL, FTL_REFRESH_PERIOD) *
                     1000000000ULL);*/
      engine.scheduleEvent(
          refreshEvent,
          tick + 1800000000000000);
    });
    engine.scheduleEvent(
        refreshEvent, 1800000000000000);
  }

  // Refresh static file
  //refreshStatFile.open("/home/rdolf/EE817/simplessd/log/refresh_stat_10^-3_1hour_no_accel.txt");
  stat.refreshCallCount = 0;
  debugprint(LOG_FTL_PAGE_MAPPING, "Refresh setting done. The number of queues: %u", refreshQueues.size());

  // Step 1. Filling
  if (mode == FILLING_MODE_0 || mode == FILLING_MODE_1) {
    // Sequential
    for (uint64_t i = 0; i < nPagesToWarmup; i++) {
      tick = 0;
      req.lpn = i;
      writeInternal(req, tick, false);
    }
  }
  else {
    // Random
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(0, nTotalLogicalPages - 1);

    for (uint64_t i = 0; i < nPagesToWarmup; i++) {
      tick = 0;
      req.lpn = dist(gen);
      writeInternal(req, tick, false);
    }
  }

  // Step 2. Invalidating
  if (mode == FILLING_MODE_0) {
    // Sequential
    for (uint64_t i = 0; i < nPagesToInvalidate; i++) {
      tick = 0;
      req.lpn = i;
      writeInternal(req, tick, false);
    }
  }
  else if (mode == FILLING_MODE_1) {
    // Random
    // We can successfully restrict range of LPN to create exact number of
    // invalid pages because we wrote in sequential mannor in step 1.
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(0, nPagesToWarmup - 1);

    for (uint64_t i = 0; i < nPagesToInvalidate; i++) {
      tick = 0;
      req.lpn = dist(gen);
      writeInternal(req, tick, false);
    }
  }
  else {
    // Random
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(0, nTotalLogicalPages - 1);

    for (uint64_t i = 0; i < nPagesToInvalidate; i++) {
      tick = 0;
      req.lpn = dist(gen);
      writeInternal(req, tick, false);
    }
  }



  // Report
  calculateTotalPages(valid, invalid);
  debugprint(LOG_FTL_PAGE_MAPPING, "Filling finished. Page status:");
  debugprint(LOG_FTL_PAGE_MAPPING,
             "  Total valid physical pages: %" PRIu64
             " (%.2f %%, target: %" PRIu64 ", error: %" PRId64 ")",
             valid, valid * 100.f / nTotalLogicalPages, nPagesToWarmup,
             (int64_t)(valid - nPagesToWarmup));
  debugprint(LOG_FTL_PAGE_MAPPING,
             "  Total invalid physical pages: %" PRIu64
             " (%.2f %%, target: %" PRIu64 ", error: %" PRId64 ")",
             invalid, invalid * 100.f / nTotalLogicalPages, nPagesToInvalidate,
             (int64_t)(invalid - nPagesToInvalidate));
  debugprint(LOG_FTL_PAGE_MAPPING, "Initialization finished");

  return true;
}


void PageMapping::refresh_event(uint64_t tick){

  uint32_t refreshMode = 1;   // 0 : Linear queue (similar with page), 1 : Scan all blocks
  float maxRBER = conf.readFloat(CONFIG_FTL, FTL_REFRESH_MAX_RBER);
  debugprint(LOG_FTL_PAGE_MAPPING, "maxRBER: %f",  maxRBER);

  if (refreshMode == 0) {
    uint64_t refreshcallCount = stat.refreshCallCount / 2;
    uint32_t target_queue = (refreshcallCount + 1) % refreshQueues.size();
    
    debugprint(LOG_FTL_PAGE_MAPPING, "Refresh call count: %lu", refreshcallCount);
    debugprint(LOG_FTL_PAGE_MAPPING, "Refresh checkPeriod: %lu", refreshcallCount);
    
    debugprint(LOG_FTL_PAGE_MAPPING, "Target queue: %u", target_queue);

    std::deque<uint32_t> tempQueue = checkedQueues[target_queue];
    checkedQueues[target_queue] = refreshQueues[target_queue];
    refreshQueues[target_queue] = tempQueue;
    refreshBlock(target_queue, tick);

    stat.refreshCallCount++;
    stat.layerCheckCount += 0;
  }
  else if (refreshMode == 1) {
    //std::vector<uint32_t> refreshTargetQueue;
    for (auto &iter : blocks) {   // Scan all blocks
      if (errorModel.getRBER(tick - iter.second.getLastWrittenTime(), iter.second.getEraseCount(), 63)
          < maxRBER  || iter.second.getInRefreshQueue() == true) {  //< 0.00032) {
        continue;
      }
      // Refresh all blocks having data retention time exceeding thredhold
      refreshTargetQueue.push_back(iter.first);
      iter.second.setInRefreshQueue(true);
    }
    //debugprint(LOG_FTL_PAGE_MAPPING, "Refresh call count: %lu",  stat.refreshCallCount);
    doRefresh(tick);

    stat.refreshCallCount++;
  }
  else if (refreshMode == 2) {
    //std::vector<uint32_t> refreshTargetQueue;
    for (auto &iter : blocks) {   // Scan all blocks
      uint32_t refreshPeriod = conf.readUint(CONFIG_FTL, FTL_REFRESH_THRESHOLD);
      if ((tick - iter.second.getLastWrittenTime()) < 3600000000000000*24*refreshPeriod
          || iter.second.getInRefreshQueue() == true) {  //< 0.00032) {
        continue;
      }
      // Refresh all blocks having data retention time exceeding thredhold
      refreshTargetQueue.push_back(iter.first);
      iter.second.setInRefreshQueue(true);
    }
    //debugprint(LOG_FTL_PAGE_MAPPING, "Refresh call count: %lu",  stat.refreshCallCount);
    doRefresh(tick);

    stat.refreshCallCount++;
  }
}



// insert to refresh queue
void PageMapping::setRefreshPeriod(uint32_t eraseCount, uint32_t blockID, uint32_t layerID){

  uint32_t refreshMode = 1;
  if (refreshMode != 0){
    return;
  }
  

  uint64_t refreshcallCount = stat.refreshCallCount / 2;
  uint32_t num_queue = refreshQueues.size();
  uint32_t cur_queue = refreshcallCount % num_queue;
  
  auto block = blocks.find(blockID);
  //debugprint(LOG_FTL_PAGE_MAPPING, "Set refresh start: %u, %u", eraseCount, layerID);
  for (uint32_t i = 1, j = 1; i <= num_queue; i++, j=j+1){
    if (i == num_queue) {
      refreshQueues[cur_queue].push_back(blockID);
      block->second.setRefreshQueueNum(cur_queue);
      //refreshQueues[i-1].push_back(blockID);
      //block->second.setRefreshQueueNum(i-1);
      break;
    }
    //std::cout << "j " << j << std::endl;
    //debugprint(LOG_FTL_PAGE_MAPPING, "refresh period: %lu", refresh_period);
    float newRBER = errorModel.getRBER(refresh_period * j, eraseCount, layerID);
    //debugprint(LOG_FTL_PAGE_MAPPING, "%u period RBER: %f", i, newRBER);

    if (newRBER > 0.00032){ // 10^-4 = ECC capability
      //debugprint(LOG_FTL_PAGE_MAPPING, "insert %u, %u, %u", block->first, layerNumber, i);
      
      
      if ((block->second.getRefreshQueueNum() - cur_queue) > j) {
        refreshQueues[(cur_queue + j) % num_queue].push_back(blockID);
        block->second.setRefreshQueueNum((cur_queue + j) % num_queue);
        
      }
  
      break; // In deque version, we should insert layer to only one deque
    }
  }
  insertedBlockCheck.set(blockID,true);
  //debugprint(LOG_FTL_PAGE_MAPPING, "Set refresh end: %u, %u", eraseCount, layerID);
}


void PageMapping::read(Request &req, uint64_t &tick) {
  uint64_t begin = tick;

  if (req.ioFlag.count() > 0) {
    readInternal(req, tick);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "READ  | LPN %" PRIu64 " | %" PRIu64 " - %" PRIu64 " (%" PRIu64
               ")",
               req.lpn, begin, tick, tick - begin);
  }
  else {
    warn("FTL got empty request");
  }

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::READ);
}

void PageMapping::write(Request &req, uint64_t &tick) {
  uint64_t begin = tick;
  
  if (req.ioFlag.count() > 0) {
    writeInternal(req, tick);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "WRITE | LPN %" PRIu64 " | %" PRIu64 " - %" PRIu64 " (%" PRIu64
               ")",
               req.lpn, begin, tick, tick - begin);
  }
  else {
    warn("FTL got empty request");
  }

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::WRITE);
}

void PageMapping::trim(Request &req, uint64_t &tick) {
  uint64_t begin = tick;

  trimInternal(req, tick);

  debugprint(LOG_FTL_PAGE_MAPPING,
             "TRIM  | LPN %" PRIu64 " | %" PRIu64 " - %" PRIu64 " (%" PRIu64
             ")",
             req.lpn, begin, tick, tick - begin);

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::TRIM);
}

void PageMapping::format(LPNRange &range, uint64_t &tick) {
  PAL::Request req(param.ioUnitInPage);
  std::vector<uint32_t> list;

  req.ioFlag.set();

  for (auto iter = table.begin(); iter != table.end();) {
    if (iter->first >= range.slpn && iter->first < range.slpn + range.nlp) {
      auto &mappingList = iter->second;

      // Do trim
      for (uint32_t idx = 0; idx < bitsetSize; idx++) {
        auto &mapping = mappingList.at(idx);
        auto block = blocks.find(mapping.first);

        if (block == blocks.end()) {
          panic("Block is not in use");
        }

        block->second.invalidate(mapping.second, idx);

        // Collect block indices
        list.push_back(mapping.first);
      }

      iter = table.erase(iter);
    }
    else {
      iter++;
    }
  }

  // Get blocks to erase
  std::sort(list.begin(), list.end());
  auto last = std::unique(list.begin(), list.end());
  list.erase(last, list.end());

  // Do GC only in specified blocks
  doGarbageCollection(list, tick, false);

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::FORMAT);
}

Status *PageMapping::getStatus(uint64_t lpnBegin, uint64_t lpnEnd) {
  status.freePhysicalBlocks = nFreeBlocks;

  if (lpnBegin == 0 && lpnEnd >= status.totalLogicalPages) {
    status.mappedLogicalPages = table.size();
  }
  else {
    status.mappedLogicalPages = 0;

    for (uint64_t lpn = lpnBegin; lpn < lpnEnd; lpn++) {
      if (table.count(lpn) > 0) {
        status.mappedLogicalPages++;
      }
    }
  }

  return &status;
}

float PageMapping::freeBlockRatio() {
  return (float)nFreeBlocks / param.totalPhysicalBlocks;
}

uint32_t PageMapping::convertBlockIdx(uint32_t blockIdx) {
  return blockIdx % param.pageCountToMaxPerf;
}

uint32_t PageMapping::getFreeBlock(uint32_t idx) {
  uint32_t blockIndex = 0;

  if (idx >= param.pageCountToMaxPerf) {
    panic("Index out of range");
  }

  if (nFreeBlocks > 0) {
    // Search block which is blockIdx % param.pageCountToMaxPerf == idx
    auto iter = freeBlocks.begin();

    for (; iter != freeBlocks.end(); iter++) {
      blockIndex = iter->getBlockIndex();

      if (blockIndex % param.pageCountToMaxPerf == idx) {
        break;
      }
    }

    // Sanity check
    if (iter == freeBlocks.end()) {
      // Just use first one
      iter = freeBlocks.begin();
      blockIndex = iter->getBlockIndex();
    }

    // Insert found block to block list
    if (blocks.find(blockIndex) != blocks.end()) {
      panic("Corrupted");
    }

    blocks.emplace(blockIndex, std::move(*iter));

    // Update first write time
    blocks.find(blockIndex)->second.setLastWrittenTime(getTick());
    blocks.find(blockIndex)->second.setInRefreshQueue(false);


    // Remove found block from free block list
    freeBlocks.erase(iter);
    nFreeBlocks--;
  }
  else {
    std::cout << "nFreeBlock" << nFreeBlocks << std::endl;
    panic("No free block left");
  }

  return blockIndex;
}

uint32_t PageMapping::getLastFreeBlock(Bitset &iomap) {
  if (!bRandomTweak || (lastFreeBlockIOMap & iomap).any()) {
    // Update lastFreeBlockIndex
    lastFreeBlockIndex++;

    if (lastFreeBlockIndex == param.pageCountToMaxPerf) {
      lastFreeBlockIndex = 0;
    }

    lastFreeBlockIOMap = iomap;
  }
  else {
    lastFreeBlockIOMap |= iomap;
  }

  auto freeBlock = blocks.find(lastFreeBlock.at(lastFreeBlockIndex));

  // Sanity check
  if (freeBlock == blocks.end()) {
    panic("Corrupted");
  }

  // If current free block is full, get next block
  if (freeBlock->second.getNextWritePageIndex() == param.pagesInBlock) {
    lastFreeBlock.at(lastFreeBlockIndex) = getFreeBlock(lastFreeBlockIndex);

    bReclaimMore = true;
  }

  return lastFreeBlock.at(lastFreeBlockIndex);
}

// calculate weight of each block regarding victim selection policy
void PageMapping::calculateVictimWeight(
    std::vector<std::pair<uint32_t, float>> &weight, const EVICT_POLICY policy,
    uint64_t tick) {
  float temp;

  weight.reserve(blocks.size());

  switch (policy) {
    case POLICY_GREEDY:
    case POLICY_RANDOM:
    case POLICY_DCHOICE:
      for (auto &iter : blocks) {
        if (iter.second.getNextWritePageIndex() != param.pagesInBlock) {
          continue;
        }

        weight.push_back({iter.first, iter.second.getValidPageCountRaw()});
      }

      break;
    case POLICY_COST_BENEFIT:
      for (auto &iter : blocks) {
        if (iter.second.getNextWritePageIndex() != param.pagesInBlock) {
          continue;
        }

        temp = (float)(iter.second.getValidPageCountRaw()) / param.pagesInBlock;

        weight.push_back(
            {iter.first,
             temp / ((1 - temp) * (tick - iter.second.getLastAccessedTime()))});
      }

      break;
    default:
      panic("Invalid evict policy");
  }
}

void PageMapping::selectVictimBlock(std::vector<uint32_t> &list,
                                    uint64_t &tick, std::vector<uint32_t> &exceptList) {
  static const GC_MODE mode = (GC_MODE)conf.readInt(CONFIG_FTL, FTL_GC_MODE);
  static const EVICT_POLICY policy =
      (EVICT_POLICY)conf.readInt(CONFIG_FTL, FTL_GC_EVICT_POLICY);
  static uint32_t dChoiceParam =
      conf.readUint(CONFIG_FTL, FTL_GC_D_CHOICE_PARAM);
  uint64_t nBlocks = conf.readUint(CONFIG_FTL, FTL_GC_RECLAIM_BLOCK);
  std::vector<std::pair<uint32_t, float>> weight;

  list.clear();

  // Calculate number of blocks to reclaim
  if (mode == GC_MODE_0) {
    // DO NOTHING
  }
  else if (mode == GC_MODE_1) {
    static const float t = conf.readFloat(CONFIG_FTL, FTL_GC_RECLAIM_THRESHOLD);

    nBlocks = param.totalPhysicalBlocks * t - nFreeBlocks;
  }
  else {
    panic("Invalid GC mode");
  }

  // reclaim one more if last free block fully used
  if (bReclaimMore) {
    nBlocks += param.pageCountToMaxPerf;

    bReclaimMore = false;
  }

  // Calculate weights of all blocks
  calculateVictimWeight(weight, policy, tick);

  if (policy == POLICY_RANDOM || policy == POLICY_DCHOICE) {
    uint64_t randomRange =
        policy == POLICY_RANDOM ? nBlocks : dChoiceParam * nBlocks;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(0, weight.size() - 1);
    std::vector<std::pair<uint32_t, float>> selected;

    while (selected.size() < randomRange) {
      uint64_t idx = dist(gen);

      auto findIter = std::find(exceptList.begin(), exceptList.end(), weight.at(idx).first);

      if (weight.at(idx).first < std::numeric_limits<uint32_t>::max() 
          && findIter == exceptList.end()) {
        selected.push_back(weight.at(idx));
        weight.at(idx).first = std::numeric_limits<uint32_t>::max();
      }
    }

    weight = std::move(selected);
  }

  // Sort weights
  std::sort(
      weight.begin(), weight.end(),
      [](std::pair<uint32_t, float> a, std::pair<uint32_t, float> b) -> bool {
        return a.second < b.second;
      });

  // Select victims from the blocks with the lowest weight
  nBlocks = MIN(nBlocks, weight.size());

  for (uint64_t i = 0; i < nBlocks; i++) {
    list.push_back(weight.at(i).first);
  }

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::SELECT_VICTIM_BLOCK);
}

void PageMapping::doGarbageCollection(std::vector<uint32_t> &blocksToReclaim,
                                      uint64_t &tick, bool isRefresh) {
  PAL::Request req(param.ioUnitInPage);
  std::vector<PAL::Request> readRequests;
  std::vector<PAL::Request> writeRequests;
  std::vector<PAL::Request> eraseRequests;
  std::vector<uint64_t> lpns;
  Bitset bit(param.ioUnitInPage);
  uint64_t beginAt;
  uint64_t readFinishedAt = tick;
  uint64_t writeFinishedAt = tick;
  uint64_t eraseFinishedAt = tick;

  if (blocksToReclaim.size() == 0) {
    return;
  }

  // For all blocks to reclaim, collecting request structure only
  for (auto &iter : blocksToReclaim) {
    auto block = blocks.find(iter);

    if (block == blocks.end()) {
      panic("Invalid block");
    }

    // Copy valid pages to free block
    for (uint32_t pageIndex = 0; pageIndex < param.pagesInBlock; pageIndex++) {
      // Valid?
      if (block->second.getPageInfo(pageIndex, lpns, bit)) {
        if (!bRandomTweak) {
          bit.set();
        }

        // Retrive free block
        auto freeBlock = blocks.find(getLastFreeBlock(bit));

        // Issue Read
        req.blockIndex = block->first;
        req.pageIndex = pageIndex;
        req.ioFlag = bit;

        readRequests.push_back(req);

        // Update mapping table
        uint32_t newBlockIdx = freeBlock->first;

        for (uint32_t idx = 0; idx < bitsetSize; idx++) {
          if (bit.test(idx)) {
            // Invalidate
            block->second.invalidate(pageIndex, idx);

            auto mappingList = table.find(lpns.at(idx));

            if (mappingList == table.end()) {
              panic("Invalid mapping table entry");
            }

            pDRAM->read(&(*mappingList), 8 * param.ioUnitInPage, tick);

            auto &mapping = mappingList->second.at(idx);

            uint32_t newPageIdx = freeBlock->second.getNextWritePageIndex(idx);

            mapping.first = newBlockIdx;
            mapping.second = newPageIdx;

            freeBlock->second.write(newPageIdx, lpns.at(idx), idx, beginAt);

            // Issue Write
            req.blockIndex = newBlockIdx;
            req.pageIndex = newPageIdx;

            if (bRandomTweak) {
              req.ioFlag.reset();
              req.ioFlag.set(idx);
            }
            else {
              req.ioFlag.set();
            }

            writeRequests.push_back(req);

            // set new refresh period
            uint32_t eraseCount = freeBlock->second.getEraseCount();
            uint32_t layerNumber = newPageIdx % 64;

            //debugprint(LOG_FTL_PAGE_MAPPING, "set refresh period - erasecount, layerNymber, blockIdx, pageIdx: %u, %u, %u, %u",
            //          eraseCount, layerNumber, newBlockIdx, newPageIdx);
            setRefreshPeriod(eraseCount, newBlockIdx, layerNumber);

            if(isRefresh){
              stat.refreshPageCopies++;
            }
            stat.validPageCopies++;
          }
        }
        // TODO : this will updated for every write. The refresh should be done only for first write
        //freeBlock->second.setLastWrittenTime(tick);   
        stat.validSuperPageCopies++;
      }
    }

    // Erase block
    req.blockIndex = block->first;
    req.pageIndex = 0;
    req.ioFlag.set();

    eraseRequests.push_back(req);
  }

  // Do actual I/O here
  // This handles PAL2 limitation (SIGSEGV, infinite loop, or so-on)
  
  for (auto &iter : readRequests) {
    beginAt = tick;

    pPAL->read(iter, beginAt);

    readFinishedAt = MAX(readFinishedAt, beginAt);
  }

  for (auto &iter : writeRequests) {
    beginAt = readFinishedAt;

    pPAL->write(iter, beginAt);

    writeFinishedAt = MAX(writeFinishedAt, beginAt);
  }

  for (auto &iter : eraseRequests) {
    beginAt = readFinishedAt;

    eraseInternal(iter, beginAt);

    eraseFinishedAt = MAX(eraseFinishedAt, beginAt);
  }
  


  tick = MAX(writeFinishedAt, eraseFinishedAt);
  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::DO_GARBAGE_COLLECTION);
}


void PageMapping::doRefresh(uint64_t &tick) {
  PAL::Request req(param.ioUnitInPage);
  std::vector<PAL::Request> readRequests;
  std::vector<PAL::Request> writeRequests;
  std::vector<PAL::Request> eraseRequests;
  std::vector<uint64_t> lpns;
  Bitset bit(param.ioUnitInPage);
  uint64_t beginAt;
  uint64_t readFinishedAt = tick;
  uint64_t writeFinishedAt = tick;
  uint64_t eraseFinishedAt = tick;

  std::vector<uint64_t> tempLpns;
  Bitset tempBit(param.ioUnitInPage);
  static float gcThreshold = conf.readFloat(CONFIG_FTL, FTL_GC_THRESHOLD_RATIO);

  if (refreshTargetQueue.size() == 0) {
    return;
  }


  if (freeBlockRatio() < gcThreshold) {
    
    debugprint(LOG_FTL_PAGE_MAPPING, "gcThreshold : %lf", gcThreshold);
    debugprint(LOG_FTL_PAGE_MAPPING, "freeBlockRatio : %lf", freeBlockRatio());
    debugprint(LOG_FTL_PAGE_MAPPING, "n free blocks : %u", nFreeBlocks);

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    std::vector<uint32_t> dummy;

    selectVictimBlock(list, beginAt, dummy);

    // If the block would be garbage collected, it shouldn't be refeshed
    //for (auto & gcIter : list) {
      //debugprint(LOG_FTL_PAGE_MAPPING, "Block %u will be garbage collected", gcIter);
    //  blocksToRefresh.erase(std::remove(blocksToRefresh.begin(), blocksToRefresh.end(), gcIter), blocksToRefresh.end());
    //}
    

    debugprint(LOG_FTL_PAGE_MAPPING,
              "GC   | Refreshing | %u blocks will be reclaimed", list.size());

    doGarbageCollection(list, beginAt, true);

    debugprint(LOG_FTL_PAGE_MAPPING,
              "GC   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
              beginAt, beginAt - tick);
    stat.gcCount++;
    stat.reclaimedBlocks += list.size();
    //debugprint(LOG_FTL_PAGE_MAPPING, "n free blocks after gc : %u", nFreeBlocks);

    // Problem : refresh 될 block이 garbage collection에 의해 erase 될 수 있음
  }
  
  uint32_t maxRefreshBlock = 1000;
  //debugprint(LOG_FTL_PAGE_MAPPING, "start refreshing");
  // For all blocks to reclaim, collecting request structure only
  for (uint32_t i = 0; i < maxRefreshBlock; i++){

    if (refreshTargetQueue.size() == 0) {
      break;
    }
    uint32_t blockIndex = refreshTargetQueue.front();
    refreshTargetQueue.pop_front();
    
    auto block = blocks.find(blockIndex);

    if (block == blocks.end()) {
      //printf("Cannot find block %u", iter);
      //panic("Invalid block, refresh failed");
      // This can be happened when the block is GCed before refresh
      i = i -1;
      continue;
    }

    bool refreshCheck = false;

    // Copy valid pages to free block
    for (uint32_t pageIndex = 0; pageIndex < param.pagesInBlock; pageIndex++) {
      //debugprint(LOG_FTL_PAGE_MAPPING, "Check valid");
      // Valid?
      if (block->second.getPageInfo(pageIndex, lpns, bit)) {
        if (!bRandomTweak) {
          bit.set();
        }

        if (refreshCheck == false){
          refreshCheck = true;
        }
        

        // Retrive free block
        auto freeBlock = blocks.find(getLastFreeBlock(bit));

        // Issue Read
        req.blockIndex = block->first;
        req.pageIndex = pageIndex;
        req.ioFlag = bit;

        readRequests.push_back(req);

        // Update mapping table
        uint32_t newBlockIdx = freeBlock->first;

        for (uint32_t idx = 0; idx < bitsetSize; idx++) {
          if (bit.test(idx)) {    
            // Invalidate
            block->second.invalidate(pageIndex, idx); // 여기서 out of range error

            auto mappingList = table.find(lpns.at(idx));

            if (mappingList == table.end()) {
              panic("Invalid mapping table entry, refresh failed");
            }

            pDRAM->read(&(*mappingList), 8 * param.ioUnitInPage, tick);

            auto &mapping = mappingList->second.at(idx);

            uint32_t newPageIdx = freeBlock->second.getNextWritePageIndex(idx);

            mapping.first = newBlockIdx;
            mapping.second = newPageIdx;


            freeBlock->second.write(newPageIdx, lpns.at(idx), idx, beginAt);

            //freeBlock->second.getPageInfo(newPageIdx, tempLpns, tempBit);
            //debugprint(LOG_FTL_PAGE_MAPPING, "got page info");
 
            // Issue Write
            req.blockIndex = newBlockIdx;
            req.pageIndex = newPageIdx;

            if (bRandomTweak) {
              req.ioFlag.reset();
              req.ioFlag.set(idx);
            }
            else {
              req.ioFlag.set();
            }

            writeRequests.push_back(req);

            stat.refreshPageCopies++;
          }
        }
        stat.refreshSuperPageCopies++;
      }
    }
    // TODO: Should be garbage collected when there is not enough blocks
    // Or write should be performed by writeInternal
    //stat.refreshedBlocks++;
    if (refreshCheck == true){
      stat.refreshedBlocks++;
    }
  }
  //debugprint(LOG_FTL_PAGE_MAPPING, "Do actual I/O");
  // Do actual I/O here
  // This handles PAL2 limitation (SIGSEGV, infinite loop, or so-on)
  for (auto &iter : readRequests) {
    beginAt = tick;

    pPAL->read(iter, beginAt);

    readFinishedAt = MAX(readFinishedAt, beginAt);
  }

  for (auto &iter : writeRequests) {
    beginAt = readFinishedAt;

    pPAL->write(iter, beginAt);

    writeFinishedAt = MAX(writeFinishedAt, beginAt);
  }

  tick = MAX(writeFinishedAt, eraseFinishedAt);
  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::DO_GARBAGE_COLLECTION);

  if (freeBlockRatio() < gcThreshold) {
    
    debugprint(LOG_FTL_PAGE_MAPPING, "gcThreshold : %lf", gcThreshold);
    debugprint(LOG_FTL_PAGE_MAPPING, "freeBlockRatio : %lf", freeBlockRatio());
    debugprint(LOG_FTL_PAGE_MAPPING, "n free blocks : %u", nFreeBlocks);

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    std::vector<uint32_t> dummy;

    selectVictimBlock(list, beginAt, dummy);

    // If the block would be garbage collected, it shouldn't be refeshed
    //for (auto & gcIter : list) {
      //debugprint(LOG_FTL_PAGE_MAPPING, "Block %u will be garbage collected", gcIter);
    //  blocksToRefresh.erase(std::remove(blocksToRefresh.begin(), blocksToRefresh.end(), gcIter), blocksToRefresh.end());
    //}
    

    debugprint(LOG_FTL_PAGE_MAPPING,
              "GC   | Refreshing | %u blocks will be reclaimed", list.size());

    doGarbageCollection(list, beginAt, true);

    debugprint(LOG_FTL_PAGE_MAPPING,
              "GC   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
              beginAt, beginAt - tick);
    stat.gcCount++;
    stat.reclaimedBlocks += list.size();
    //debugprint(LOG_FTL_PAGE_MAPPING, "n free blocks after gc : %u", nFreeBlocks);

    // Problem : refresh 될 block이 garbage collection에 의해 erase 될 수 있음
  }
  
}

void PageMapping::refreshBlock(uint32_t queueNum, uint64_t &tick) {
  //debugprint(LOG_FTL_PAGE_MAPPING, "Refresh page start");
  PAL::Request req(param.ioUnitInPage);
  std::vector<PAL::Request> readRequests;
  std::vector<PAL::Request> writeRequests;
  std::vector<PAL::Request> eraseRequests;
  std::vector<uint64_t> lpns;
  Bitset bit(param.ioUnitInPage);
  uint64_t beginAt = tick;
  uint64_t readFinishedAt = tick;
  uint64_t writeFinishedAt = tick;
  //uint64_t eraseFinishedAt = tick;

  std::vector<uint64_t> tempLpns;
  Bitset tempBit(param.ioUnitInPage);
  static float gcThreshold = conf.readFloat(CONFIG_FTL, FTL_GC_THRESHOLD_RATIO);

  // GC before refresh
  if (freeBlockRatio() < gcThreshold) {

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    std::vector<uint32_t> dummy;
    selectVictimBlock(list, beginAt, dummy);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "GC   | Refreshing | %u blocks will be reclaimed", list.size());

    doGarbageCollection(list, beginAt, true);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "GC   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
               beginAt, beginAt - tick);

    stat.gcCount++;
    stat.reclaimedBlocks += list.size();
  }
  
  // For all blocks to reclaim, collecting request structure only
  uint32_t maxRefreshBlock = 6000;  // # of blocks can be refreshed for each refresh interval (not check interval)

  for (uint32_t i = 0; i < maxRefreshBlock; i++) {
    if (checkedQueues[queueNum].empty()){
      break;
    }
    
    uint32_t blockIndex = checkedQueues[queueNum].front();
    checkedQueues[queueNum].pop_front();

    if (!insertedBlockCheck.test(blockIndex)) {
      // This can be happened, when block is refreshed already with high priority
      // Block can be inserted multiple times e.g. inserted by layer with low RBER, then inserted by layer with high RBER
      i = i -1;
      continue;
      //panic("Corrupted block check bitmap, refresh failed. Block ID : %u", blockIndex);
    }

    // Reset refresh check bitmap
    insertedBlockCheck.set(blockIndex,false);
    
    auto block = blocks.find(blockIndex);
    if (block == blocks.end()) {
      //panic("Invalid block, refresh failed");
      // This can be happen if the block is GCed at the beginning of refresh
      i = i -1;
      continue;
    }

      // Copy valid pages to free block
      for (uint32_t pageIndex = 0; pageIndex < param.pagesInBlock; pageIndex ++) {

        //if (block->second.getValidPageCount()) {  // Valid?
        if (block->second.getPageInfo(pageIndex, lpns, bit)) {  //Modified!!!
          if (!bRandomTweak) {
            bit.set();
          }

          // Retrive free block
          auto freeBlock = blocks.find(getLastFreeBlock(bit));

          // Issue Read
          req.blockIndex = block->first;
          req.pageIndex = pageIndex;
          req.ioFlag = bit;

          readRequests.push_back(req);

          // Update mapping table
          uint32_t newBlockIdx = freeBlock->first;

          for (uint32_t idx = 0; idx < bitsetSize; idx++) {
            if (bit.test(idx)) {    
              // Invalidate
              block->second.invalidate(pageIndex, idx); // 여기서 out of range error 났었음 (왜났었는지는 기억이 안남)

              auto mappingList = table.find(lpns.at(idx));

              if (mappingList == table.end()) {     //Modified!!!
                panic("Invalid mapping table entry, refresh failed");
                // This shouldn't be happened. For all valid page, there should be valid mapping
                //continue;
              }

              pDRAM->read(&(*mappingList), 8 * param.ioUnitInPage, tick);

              auto &mapping = mappingList->second.at(idx);

              uint32_t newPageIdx = freeBlock->second.getNextWritePageIndex(idx);

              mapping.first = newBlockIdx;
              mapping.second = newPageIdx;

              freeBlock->second.write(newPageIdx, lpns.at(idx), idx, beginAt);

              //freeBlock->second.getPageInfo(newPageIdx, tempLpns, tempBit);
              //debugprint(LOG_FTL_PAGE_MAPPING, "got page info");

              // Issue Write
              req.blockIndex = newBlockIdx;
              req.pageIndex = newPageIdx;

              if (bRandomTweak) {
                req.ioFlag.reset();
                req.ioFlag.set(idx);
              }
              else {
                req.ioFlag.set();
              }

              writeRequests.push_back(req);
              
              uint32_t eraseCount = freeBlock->second.getEraseCount();
              uint32_t layerNumber = newPageIdx % 64;

              setRefreshPeriod(eraseCount, newBlockIdx, layerNumber);

              stat.refreshPageCopies++;
            }
          }
          //debugprint(LOG_FTL_PAGE_MAPPING, "set last written time");
          //freeBlock->second.setLastWrittenTime(tick);

          stat.refreshSuperPageCopies++;
        }
      }
    stat.refreshedBlocks++;
  }


  // Copy valid pages to free block
  // pageIndex is n*layerNum??
  //layerNum = 0;
  
  // TODO: Should be garbage collected when there is not enough blocks
  // Or write should be performed by writeInternal

  //debugprint(LOG_FTL_PAGE_MAPPING, "Do actual I/O");
  // Do actual I/O here
  // This handles PAL2 limitation (SIGSEGV, infinite loop, or so-on)
  for (auto &iter : readRequests) {
    beginAt = tick;

    pPAL->read(iter, beginAt);

    readFinishedAt = MAX(readFinishedAt, beginAt);
  }

  for (auto &iter : writeRequests) {
    beginAt = readFinishedAt;

    pPAL->write(iter, beginAt);

    writeFinishedAt = MAX(writeFinishedAt, beginAt);
  }
  
  //debugprint(LOG_FTL_PAGE_MAPPING, "page refresh done. remaining free blocks: %u", nFreeBlocks);
  tick = MAX(writeFinishedAt, readFinishedAt);
  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::DO_GARBAGE_COLLECTION);
    // GC after refresh
  if (freeBlockRatio() < gcThreshold) {

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    std::vector<uint32_t> dummy;
    selectVictimBlock(list, beginAt, dummy);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "GC   | Refreshing | %u blocks will be reclaimed", list.size());

    doGarbageCollection(list, beginAt, true);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "GC   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
               beginAt, beginAt - tick);

    stat.gcCount++;
    stat.reclaimedBlocks += list.size();
  }
  //debugprint(LOG_FTL_PAGE_MAPPING, "Refresh page end");
}


void PageMapping::readInternal(Request &req, uint64_t &tick) {
  PAL::Request palRequest(req);
  uint64_t beginAt;
  uint64_t finishedAt = tick;

  auto mappingList = table.find(req.lpn);

  if (mappingList != table.end()) {
    if (bRandomTweak) {
      pDRAM->read(&(*mappingList), 8 * req.ioFlag.count(), tick);
    }
    else {
      pDRAM->read(&(*mappingList), 8, tick);
    }

    for (uint32_t idx = 0; idx < bitsetSize; idx++) {
      if (req.ioFlag.test(idx) || !bRandomTweak) {
        auto &mapping = mappingList->second.at(idx);

        if (mapping.first < param.totalPhysicalBlocks &&
            mapping.second < param.pagesInBlock) {
          palRequest.blockIndex = mapping.first;
          palRequest.pageIndex = mapping.second;

          if (bRandomTweak) {
            palRequest.ioFlag.reset();
            palRequest.ioFlag.set(idx);
          }
          else {
            palRequest.ioFlag.set();
          }

          auto block = blocks.find(palRequest.blockIndex);

          if (block == blocks.end()) {
            panic("Block is not in use");
          }

          beginAt = tick;

          block->second.read(palRequest.pageIndex, idx, beginAt);
          pPAL->read(palRequest, beginAt);

          /*
          uint64_t lastWritten = block->second.getLastWrittenTime();
          uint32_t eraseCount = block->second.getEraseCount();
          uint64_t curErrorCount = block->second.getMaxErrorCount();

          debugprint(LOG_FTL_PAGE_MAPPING, "Erase count %u", eraseCount);

          //TODO: Get layer number
          uint32_t layerNumber = mapping.second % 64;
          uint64_t newErrorCount = errorModel.getRandError(tick - lastWritten, eraseCount, layerNumber);

          debugprint(LOG_FTL_PAGE_MAPPING, "new rber: %f", errorModel.getRBER(tick - lastWritten, eraseCount, 0));
          debugprint(LOG_FTL_PAGE_MAPPING, "new randerror: %u", newErrorCount);


          block->second.setMaxErrorCount(max(curErrorCount, newErrorCount));
          */

          finishedAt = MAX(finishedAt, beginAt);
        }
      }
    }



    tick = finishedAt;
    tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::READ_INTERNAL);
  }
}

void PageMapping::writeInternal(Request &req, uint64_t &tick, bool sendToPAL) {
  //debugprint(LOG_FTL_PAGE_MAPPING, "Write internal start");
  PAL::Request palRequest(req);
  std::unordered_map<uint32_t, Block>::iterator block;
  auto mappingList = table.find(req.lpn);
  uint64_t beginAt;
  uint64_t finishedAt = tick;
  bool readBeforeWrite = false;


  if (mappingList != table.end()) {
    for (uint32_t idx = 0; idx < bitsetSize; idx++) {
      if (req.ioFlag.test(idx) || !bRandomTweak) {
        auto &mapping = mappingList->second.at(idx);

        if (mapping.first < param.totalPhysicalBlocks &&
            mapping.second < param.pagesInBlock) {
          block = blocks.find(mapping.first);

          // Invalidate current page
          block->second.invalidate(mapping.second, idx);
        }
      }
    }
  }
  else {
    // Create empty mapping
    auto ret = table.emplace(
        req.lpn,
        std::vector<std::pair<uint32_t, uint32_t>>(
            bitsetSize, {param.totalPhysicalBlocks, param.pagesInBlock}));

    if (!ret.second) {
      panic("Failed to insert new mapping");
    }

    mappingList = ret.first;
  }

  // Write data to free block
  block = blocks.find(getLastFreeBlock(req.ioFlag));

  if (block == blocks.end()) {
    panic("No such block");
  }

  if (sendToPAL) {
    if (bRandomTweak) {
      pDRAM->read(&(*mappingList), 8 * req.ioFlag.count(), tick);
      pDRAM->write(&(*mappingList), 8 * req.ioFlag.count(), tick);
    }
    else {
      pDRAM->read(&(*mappingList), 8, tick);
      pDRAM->write(&(*mappingList), 8, tick);
    }
  }

  if (!bRandomTweak && !req.ioFlag.all()) {
    // We have to read old data
    readBeforeWrite = true;
  }

  for (uint32_t idx = 0; idx < bitsetSize; idx++) {
    if (req.ioFlag.test(idx) || !bRandomTweak) {
      uint32_t pageIndex = block->second.getNextWritePageIndex(idx);
      auto &mapping = mappingList->second.at(idx);

      beginAt = tick;

      block->second.write(pageIndex, req.lpn, idx, beginAt);

      // Read old data if needed (Only executed when bRandomTweak = false)
      // Maybe some other init procedures want to perform 'partial-write'
      // So check sendToPAL variable
      if (readBeforeWrite && sendToPAL) {
        palRequest.blockIndex = mapping.first;
        palRequest.pageIndex = mapping.second;

        // We don't need to read old data
        palRequest.ioFlag = req.ioFlag;
        palRequest.ioFlag.flip();

        pPAL->read(palRequest, beginAt);
      }

      // update mapping to table
      mapping.first = block->first;
      mapping.second = pageIndex;

      if (sendToPAL) {
        palRequest.blockIndex = block->first;
        palRequest.pageIndex = pageIndex;

        if (bRandomTweak) {
          palRequest.ioFlag.reset();
          palRequest.ioFlag.set(idx);
        }
        else {
          palRequest.ioFlag.set();
        }

        pPAL->write(palRequest, beginAt);
      }

      finishedAt = MAX(finishedAt, beginAt);

      //if (sendToPAL){
      // Predict error
      uint32_t eraseCount = block->second.getEraseCount();
      uint32_t layerNumber = mapping.second % 64;
      
      //uint32_t globalLayerNum = (block->first * 64) + layerNumber;
      //debugprint(LOG_FTL_PAGE_MAPPING, "set refresh period - erasecount, globalLayerNum, blockIdx, pageIdx: %u, %u, %u, %u",
      //  eraseCount, globalLayerNum, block->first, mapping.second);
      setRefreshPeriod(eraseCount, block->first, layerNumber);
      //}

      //TODO: Now error count can be used to put layer to bloom filter

    }
  }

  // TODO: Have to record for each page
  //block->second.setLastWrittenTime(tick);

  // Exclude CPU operation when initializing
  if (sendToPAL) {
    tick = finishedAt;
    tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::WRITE_INTERNAL);
  }

  // GC if needed
  // I assumed that init procedure never invokes GC
  static float gcThreshold = conf.readFloat(CONFIG_FTL, FTL_GC_THRESHOLD_RATIO);

  if (freeBlockRatio() < gcThreshold) {
    if (!sendToPAL) {
      panic("ftl: GC triggered while in initialization");
    }

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    std::vector<uint32_t> dummy;
    selectVictimBlock(list, beginAt, dummy);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "GC   | On-demand | %u blocks will be reclaimed", list.size());

    doGarbageCollection(list, beginAt, false);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "GC   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
               beginAt, beginAt - tick);

    stat.gcCount++;
    stat.reclaimedBlocks += list.size();
  }
  /*
  if (tick - lastRefreshed > 1000000000 && sendToPAL){  // check every 1ms

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    selectRefreshVictim(list, beginAt);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "REFRESH   | On-time | %u blocks will be refreshed", list.size());

    doRefresh(list, beginAt);

    debugprint(LOG_FTL_PAGE_MAPPING,
               "REFRESH   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
               beginAt, beginAt - tick);

    if (list.size() > 0) {
      stat.refreshCount++;
      stat.refreshedBlocks += list.size();
    }
    lastRefreshed = tick;
  }
  */
 //debugprint(LOG_FTL_PAGE_MAPPING, "Write internal end");
}

void PageMapping::trimInternal(Request &req, uint64_t &tick) {
  auto mappingList = table.find(req.lpn);

  if (mappingList != table.end()) {
    if (bRandomTweak) {
      pDRAM->read(&(*mappingList), 8 * req.ioFlag.count(), tick);
    }
    else {
      pDRAM->read(&(*mappingList), 8, tick);
    }

    // Do trim
    for (uint32_t idx = 0; idx < bitsetSize; idx++) {
      auto &mapping = mappingList->second.at(idx);
      auto block = blocks.find(mapping.first);

      if (block == blocks.end()) {
        panic("Block is not in use");
      }

      block->second.invalidate(mapping.second, idx);
    }

    // Remove mapping
    table.erase(mappingList);

    tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::TRIM_INTERNAL);
  }
}

void PageMapping::eraseInternal(PAL::Request &req, uint64_t &tick) {
  static uint64_t threshold =
      conf.readUint(CONFIG_FTL, FTL_BAD_BLOCK_THRESHOLD);
  auto block = blocks.find(req.blockIndex);

  // Sanity checks
  if (block == blocks.end()) {
    panic("No such block");
  }

  if (block->second.getValidPageCount() != 0) {
    panic("There are valid pages in victim block");
  }

  // Erase block
  block->second.erase();

  pPAL->erase(req, tick);

  // Check erase count
  uint32_t erasedCount = block->second.getEraseCount();

  if (erasedCount < threshold) {
    // Reverse search
    auto iter = freeBlocks.end();

    while (true) {
      iter--;

      if (iter->getEraseCount() <= erasedCount) {
        // emplace: insert before pos
        iter++;

        break;
      }

      if (iter == freeBlocks.begin()) {
        break;
      }
    }

    // Insert block to free block list
    freeBlocks.emplace(iter, std::move(block->second));
    nFreeBlocks++;
  }

  // Remove block from block list
  blocks.erase(block);

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::ERASE_INTERNAL);
}

float PageMapping::calculateWearLeveling() {
  uint64_t totalEraseCnt = 0;
  uint64_t sumOfSquaredEraseCnt = 0;
  uint64_t numOfBlocks = param.totalLogicalBlocks;
  uint64_t eraseCnt;

  for (auto &iter : blocks) {
    eraseCnt = iter.second.getEraseCount();
    totalEraseCnt += eraseCnt;
    sumOfSquaredEraseCnt += eraseCnt * eraseCnt;
  }

  // freeBlocks is sorted
  // Calculate from backward, stop when eraseCnt is zero
  for (auto riter = freeBlocks.rbegin(); riter != freeBlocks.rend(); riter++) {
    eraseCnt = riter->getEraseCount();

    if (eraseCnt == 0) {
      break;
    }

    totalEraseCnt += eraseCnt;
    sumOfSquaredEraseCnt += eraseCnt * eraseCnt;
  }

  if (sumOfSquaredEraseCnt == 0) {
    return -1;  // no meaning of wear-leveling
  }

  return (float)totalEraseCnt * totalEraseCnt /
         (numOfBlocks * sumOfSquaredEraseCnt);
}

void PageMapping::calculateTotalPages(uint64_t &valid, uint64_t &invalid) {
  valid = 0;
  invalid = 0;

  for (auto &iter : blocks) {
    valid += iter.second.getValidPageCount();
    invalid += iter.second.getDirtyPageCount();
  }
}

float PageMapping::calculateAverageError(){
  uint64_t totalError = 0;
  float validBlockCount = 0;

  for (auto &iter : blocks) {
    totalError = totalError + iter.second.getMaxErrorCount();  
    validBlockCount = validBlockCount + 1;
  }

  float averageError = totalError / validBlockCount;

  return averageError;
}

void PageMapping::getStatList(std::vector<Stats> &list, std::string prefix) {
  //debugprint(LOG_FTL_PAGE_MAPPING, "get stat list start");
  Stats temp;

  temp.name = prefix + "page_mapping.gc.count";
  temp.desc = "Total GC count";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.gc.reclaimed_blocks";
  temp.desc = "Total reclaimed blocks in GC";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.gc.superpage_copies";
  temp.desc = "Total copied valid superpages during GC";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.gc.page_copies";
  temp.desc = "Total copied valid pages during GC";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refreshGc.page_copies";
  temp.desc = "Total copied valid pages during GC due to refresh";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refresh.count";
  temp.desc = "Total Refresh count";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refresh.refreshed_blocks";
  temp.desc = "Total blocks been refreshed";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refresh.superpage_copies";
  temp.desc = "Total copied valid superpages during Refresh";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refresh.page_copies";
  temp.desc = "Total copied valid pages during Refresh";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refresh.call_count";
  temp.desc = "The number of refresh call";
  list.push_back(temp);

  /*
  temp.name = prefix + "page_mapping.refresh.layer_check_count";
  temp.desc = "The number of total layer check";
  list.push_back(temp);
  */

  temp.name = prefix + "page_mapping.refresh.error_counts";
  temp.desc = "The average number of errors";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.refresh.refreshQueueSize";
  temp.desc = "Refresh queue size in refresh mode 1";
  list.push_back(temp);

  // For the exact definition, see following paper:
  // Li, Yongkun, Patrick PC Lee, and John Lui.
  // "Stochastic modeling of large-scale solid-state storage systems: analysis,
  // design tradeoffs and optimization." ACM SIGMETRICS (2013)
  temp.name = prefix + "page_mapping.wear_leveling";
  temp.desc = "Wear-leveling factor";
  list.push_back(temp);

  temp.name = prefix + "page_mapping.freeBlock_counts";
  temp.desc = "The number of free blocks left";
  list.push_back(temp);

  uint32_t refreshMode = 1;
  if (refreshMode == 0) {
    if (refreshQueues.size()){
      for(uint32_t i=0; i<refreshQueues.size(); i++){
        temp.name = prefix + "page_mapping.bloomFilter";
        temp.desc = "The number elements of bf-";
        list.push_back(temp);
      }
    }
  }
  

  //debugprint(LOG_FTL_PAGE_MAPPING, "get stat list end");
}

void PageMapping::getStatValues(std::vector<double> &values) {
  //debugprint(LOG_FTL_PAGE_MAPPING, "get stat values start");
  values.push_back(stat.gcCount);
  values.push_back(stat.reclaimedBlocks);
  values.push_back(stat.validSuperPageCopies);
  values.push_back(stat.validPageCopies);
  values.push_back(stat.refreshGCPageCopies);
 
  values.push_back(stat.refreshCount);
  values.push_back(stat.refreshedBlocks);
  values.push_back(stat.refreshSuperPageCopies);
  values.push_back(stat.refreshPageCopies);
  values.push_back(stat.refreshCallCount);
  //values.push_back(stat.layerCheckCount);
  values.push_back(calculateAverageError());
  values.push_back(refreshTargetQueue.size());
  values.push_back(calculateWearLeveling());

  values.push_back(nFreeBlocks);

  uint32_t refreshMode = 1;
  if (refreshMode == 0) {
    if (refreshQueues.size()){
      for(uint32_t i=0; i<refreshQueues.size(); i++){
        values.push_back(refreshQueues[i].size());
      }
    }
  }
  
  
  //debugprint(LOG_FTL_PAGE_MAPPING, "get stat values end");
}

void PageMapping::resetStatValues() {
  memset(&stat, 0, sizeof(stat));
}

// Old refresh functions
/*
void PageMapping::doRefresh(std::vector<uint32_t> &blocksToRefresh,
                                      uint64_t &tick) {
  PAL::Request req(param.ioUnitInPage);
  std::vector<PAL::Request> readRequests;
  std::vector<PAL::Request> writeRequests;
  std::vector<PAL::Request> eraseRequests;
  std::vector<uint64_t> lpns;
  Bitset bit(param.ioUnitInPage);
  uint64_t beginAt;
  uint64_t readFinishedAt = tick;
  uint64_t writeFinishedAt = tick;
  uint64_t eraseFinishedAt = tick;

  std::vector<uint64_t> tempLpns;
  Bitset tempBit(param.ioUnitInPage);
  static float gcThreshold = conf.readFloat(CONFIG_FTL, FTL_GC_THRESHOLD_RATIO);

  if (blocksToRefresh.size() == 0) {
    return;
  }


  while (nFreeBlocks < blocksToRefresh.size() * 1.5) {
    
    debugprint(LOG_FTL_PAGE_MAPPING, "gcThreshold : %lf", gcThreshold);
    debugprint(LOG_FTL_PAGE_MAPPING, "freeBlockRatio : %lf", freeBlockRatio());
    debugprint(LOG_FTL_PAGE_MAPPING, "n free blocks : %u", nFreeBlocks);

    std::vector<uint32_t> list;
    uint64_t beginAt = tick;

    std::vector<uint32_t> dummy;

    selectVictimBlock(list, beginAt, dummy);

    // If the block would be garbage collected, it shouldn't be refeshed
    for (auto & gcIter : list) {
      //debugprint(LOG_FTL_PAGE_MAPPING, "Block %u will be garbage collected", gcIter);
      blocksToRefresh.erase(std::remove(blocksToRefresh.begin(), blocksToRefresh.end(), gcIter), blocksToRefresh.end());
    }
    

    debugprint(LOG_FTL_PAGE_MAPPING,
              "GC   | Refreshing | %u blocks will be reclaimed", list.size());

    doGarbageCollection(list, beginAt);

    debugprint(LOG_FTL_PAGE_MAPPING,
              "GC   | Done | %" PRIu64 " - %" PRIu64 " (%" PRIu64 ")", tick,
              beginAt, beginAt - tick);
    stat.gcCount++;
    stat.reclaimedBlocks += list.size();
    //debugprint(LOG_FTL_PAGE_MAPPING, "n free blocks after gc : %u", nFreeBlocks);

    // Problem : refresh 될 block이 garbage collection에 의해 erase 될 수 있음
  }
  
  //debugprint(LOG_FTL_PAGE_MAPPING, "start refreshing");
  // For all blocks to reclaim, collecting request structure only
  for (auto &iter : blocksToRefresh) {
    auto block = blocks.find(iter);

    if (block == blocks.end()) {
      printf("Cannot find block %u", iter);
      panic("Invalid block, refresh failed");
    }

    // Copy valid pages to free block
    for (uint32_t pageIndex = 0; pageIndex < param.pagesInBlock; pageIndex++) {
      //debugprint(LOG_FTL_PAGE_MAPPING, "Check valid");
      // Valid?
      if (block->second.getValidPageCount()) {
        block->second.getPageInfo(pageIndex, lpns, bit);
        if (!bRandomTweak) {
          bit.set();
        }

        //debugprint(LOG_FTL_PAGE_MAPPING, "Retrive free block");

        // Retrive free block
        auto freeBlock = blocks.find(getLastFreeBlock(bit));

        // Issue Read
        req.blockIndex = block->first;
        req.pageIndex = pageIndex;
        req.ioFlag = bit;

        readRequests.push_back(req);

        //debugprint(LOG_FTL_PAGE_MAPPING, "Update mapping table");
        // Update mapping table
        uint32_t newBlockIdx = freeBlock->first;

        for (uint32_t idx = 0; idx < bitsetSize; idx++) {
          if (bit.test(idx)) {    
            //debugprint(LOG_FTL_PAGE_MAPPING, "in the if statement");
            // Invalidate
            block->second.invalidate(pageIndex, idx); // 여기서 out of range error
            //debugprint(LOG_FTL_PAGE_MAPPING, "Invalidated");


            auto mappingList = table.find(lpns.at(idx));
            //debugprint(LOG_FTL_PAGE_MAPPING, "Found mapping list");

            if (mappingList == table.end()) {
              panic("Invalid mapping table entry, refresh failed");
            }

            pDRAM->read(&(*mappingList), 8 * param.ioUnitInPage, tick);

            auto &mapping = mappingList->second.at(idx);
            //debugprint(LOG_FTL_PAGE_MAPPING, "Found mapping");

            uint32_t newPageIdx = freeBlock->second.getNextWritePageIndex(idx);

            mapping.first = newBlockIdx;
            mapping.second = newPageIdx;


            freeBlock->second.write(newPageIdx, lpns.at(idx), idx, beginAt);
            //debugprint(LOG_FTL_PAGE_MAPPING, "Written block");

            freeBlock->second.getPageInfo(newPageIdx, tempLpns, tempBit);
            //debugprint(LOG_FTL_PAGE_MAPPING, "got page info");
 
            // Issue Write
            req.blockIndex = newBlockIdx;
            req.pageIndex = newPageIdx;

            if (bRandomTweak) {
              req.ioFlag.reset();
              req.ioFlag.set(idx);
            }
            else {
              req.ioFlag.set();
            }

            writeRequests.push_back(req);

            stat.refreshPageCopies++;
          }
        }
        //debugprint(LOG_FTL_PAGE_MAPPING, "set last written time");
        //freeBlock->second.setLastWrittenTime(tick);

        stat.refreshSuperPageCopies++;
      }
    }
    // TODO: Should be garbage collected when there is not enough blocks
    // Or write should be performed by writeInternal
  }
  //debugprint(LOG_FTL_PAGE_MAPPING, "Do actual I/O");
  // Do actual I/O here
  // This handles PAL2 limitation (SIGSEGV, infinite loop, or so-on)
  for (auto &iter : readRequests) {
    beginAt = tick;

    pPAL->read(iter, beginAt);

    readFinishedAt = MAX(readFinishedAt, beginAt);
  }

  for (auto &iter : writeRequests) {
    beginAt = readFinishedAt;

    pPAL->write(iter, beginAt);

    writeFinishedAt = MAX(writeFinishedAt, beginAt);
  }

  tick = MAX(writeFinishedAt, eraseFinishedAt);
  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::DO_GARBAGE_COLLECTION);
}

// calculate weight of each block regarding victim selection policy
void PageMapping::calculateRefreshWeight(
    std::vector<std::pair<uint32_t, float>> &weight, const REFRESH_POLICY policy,
    uint64_t tick) {

  static uint64_t refreshThreshold =
      conf.readUint(CONFIG_FTL, FTL_REFRESH_THRESHOLD);

  weight.reserve(blocks.size());

  switch (policy) {
    case POLICY_NONE:
      for (auto &iter : blocks) {
        if (tick - iter.second.getLastWrittenTime() < refreshThreshold) {
          continue;
        }
        // Refresh all blocks having data retention time exceeding thredhold
        weight.push_back({iter.first, iter.second.getValidPageCountRaw()});
      }

      break;
    default:
      panic("Invalid refresh policy");
  }
}

void PageMapping::selectRefreshVictim(std::vector<uint32_t> &list,
                                    uint64_t &tick) {
  static const REFRESH_POLICY policy =
      (REFRESH_POLICY)conf.readInt(CONFIG_FTL, FTL_REFRESH_POLICY);
  std::vector<std::pair<uint32_t, float>> weight;

  list.clear();

  // Calculate weights of all blocks
  calculateRefreshWeight(weight, policy, tick);

  for (uint64_t i = 0; i < weight.size(); i++) {
    list.push_back(weight.at(i).first);
  }

  tick += applyLatency(CPU::FTL__PAGE_MAPPING, CPU::SELECT_VICTIM_BLOCK);
}
*/
}  // namespace FTL

}  // namespace SimpleSSD
