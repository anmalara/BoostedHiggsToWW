#include <unordered_map>
#include <map>
#include <vector>


template<typename TVec, typename TMap>
std::vector<TVec> RetrieveKeys(std::map<TVec, TMap> map_) {
    std::vector<TVec> vec_;
    vec_.reserve(map_.size());
    for(auto const& imap: map_) vec_.push_back(imap.first);
    return vec_;
}
