// ==========================================================================
//    SLIMM - Species Level Identification of Microbes from Metagenomes.
// ==========================================================================
// Copyright (c) 2014-2017, Temesgen H. Dadi, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#ifndef SLIMM_H
#define SLIMM_H


using namespace seqan;

inline uint32_t findLCATaxaID(std::set<uint32_t> const & taxon_ids, TNodes const & nodes);
// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class arg_options
// ----------------------------------------------------------------------------
struct arg_options
{
    typedef std::vector<std::string>            TList;

    TList rankList = {"strains",
                      "species",
                      "genus",
                      "family",
                      "order",
                      "class",
                      "phylum",
                      "superkingdom"};

    float               cov_cut_off;
    float               abundance_cut_off;
    uint32_t            bin_width;
    uint32_t            min_reads;
    bool                verbose;
    bool                is_directory;
    bool                raw_output;
    bool                coverage_output;
    std::string         rank;
    std::string         input_path;
    std::string         output_prefix;
    std::string         database_path;

    arg_options() : cov_cut_off(0.95),
                    abundance_cut_off(0.01),
                    bin_width(0),
                    min_reads(0),
                    verbose(false),
                    is_directory(false),
                    raw_output(false),
                    coverage_output(false),
                    rank("species"),
                    input_path(""),
                    output_prefix(""),
                    database_path("") {}
};

// ----------------------------------------------------------------------------
// Class slimm
// ----------------------------------------------------------------------------
class slimm
{
public:
    //constructor with argument options
    slimm(arg_options op): options(op)
    {
        collect_bam_files();
        get_considered_ranks();
        load_slimm_database(db, options.database_path);
    }

    arg_options                                         options;

    uint32_t                    current_file_index        = 0;
    uint32_t                    number_of_files           = 0;
    uint32_t                    avg_read_length           = 0;
    uint32_t                    matched_ref_length        = 0;
    uint32_t                    reference_count           = 0;
    uint32_t                    failed_by_min_read        = 0;
    uint32_t                    failed_by_min_uniq_read   = 0;
    uint32_t                    failed_byCov              = 0;
    uint32_t                    failed_byUniqCov          = 0;
    uint32_t                    hits_count                = 0;
    uint32_t                    uniq_hits_count           = 0;
    uint32_t                    matches_count             = 0;
    uint32_t                    uniq_matches_count        = 0;
    uint32_t                    uniq_matches_count2       = 0;


    slimm_database                                      db;
    std::set<uint32_t>                                  valid_ref_ids;
    std::vector<taxa_ranks>                             considered_ranks;
    std::vector<reference_contig>                       references;
    std::unordered_map<std::string, read_stat>          reads;
    std::unordered_map<uint32_t, uint32_t>              taxon_id__read_count;
    std::unordered_map<uint32_t, std::set<uint32_t> >   taxon_id__children;

    inline std::string current_bam_file_path()
    {
        return _input_paths[current_file_index];
    }

    inline void     analyze_alignments(BamFileIn & bam_file);
    inline float    coverage_cut_off();
    inline float    expected_coverage() const;
    inline void     filter_alignments();
    inline void     get_profiles();
    inline void     get_reads_lca_count();
    inline uint32_t min_reads();
    inline uint32_t min_uniq_reads();
    inline void     print_filter_stat();
    inline void     print_matches_stat();
    inline float    uniq_coverage_cut_off();
    inline void     write_raw_stat();
    inline void     write_coverage();
    inline void     write_abundance();
    inline void     reset();
    inline uint32_t get_lca(std::set<uint32_t> const & ref_ids);
    inline std::string get_lineage_string(taxa_ranks rank, std::vector<uint32_t> const & linage);
    inline std::string get_lineage_string(taxa_ranks rank, uint32_t const & taxa_id);

private:

    float                       _coverage_cut_off       = 0.0;
    float                       _uniq_coverage_cut_off  = 0.0;
    int32_t                     _min_uniq_reads         = -1;
    int32_t                     _min_reads              = -1;
    std::vector<std::string>    _input_paths;

    // member functions
    inline void collect_bam_files();
    inline void get_considered_ranks();
    inline void load_taxonomic_info();
};

inline void slimm::reset()
{
    avg_read_length           = 0;
    matched_ref_length        = 0;
    reference_count           = 0;
    failed_by_min_read        = 0;
    failed_by_min_uniq_read   = 0;
    failed_byCov              = 0;
    failed_byUniqCov          = 0;
    hits_count                = 0;
    uniq_hits_count           = 0;
    matches_count             = 0;
    uniq_matches_count        = 0;
    uniq_matches_count2       = 0;

    valid_ref_ids.clear();
    references.clear();
    reads.clear();
    taxon_id__read_count.clear();
    taxon_id__children.clear();

}


inline void slimm::analyze_alignments(BamFileIn & bam_file)
{
    BamAlignmentRecord record;
    while (!atEnd(bam_file))
    {
        readRecord(record, bam_file);
        if (hasFlagUnmapped(record) || record.rID == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.

        uint32_t center_position =  std::min(record.beginPos + (avg_read_length/2), references[record.rID].length);
        uint32_t relative_bin_no = center_position/options.bin_width;

        // maintain read properties under slimm.reads
        std::string read_name = toCString(record.qName);
        if(hasFlagFirst(record))
            append(read_name, ".1");
        else if(hasFlagLast(record))
            append(read_name, ".2");

        // if there is no read with read_name this will create one.
        reads[read_name].add_target(record.rID, relative_bin_no);
        ++hits_count;
    }


    if (hits_count == 0)
        return;

    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        if(it->second.is_uniq())
        {
            uint32_t reference_id = it->second.targets[0].reference_id;
            it->second.refs_length_sum += references[reference_id].length;
            ++uniq_matches_count;

            size_t pos_count = (it->second.targets[0]).positions.size();
            references[reference_id].reads_count += pos_count;
            it->second.refs_length_sum += references[reference_id].length;
            for (size_t j=0; j < pos_count; ++j)
            {
                uint32_t bin_number = (it->second.targets[0]).positions[j];
                ++references[reference_id].cov.bins_height[bin_number];
            }
            references[reference_id].uniq_reads_count += 1;
            uniq_hits_count += 1;
            ++references[reference_id].uniq_cov.bins_height[(it->second.targets[0]).positions[0]];
        }
        else
        {
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {

                uint32_t reference_id = it->second.targets[i].reference_id;
                it->second.refs_length_sum += references[reference_id].length;

                // ***** all of the matches in multiple pos will be counted *****
                references[reference_id].reads_count += (it->second.targets[i]).positions.size();
                for (auto bin_number : (it->second.targets[i]).positions)
                {
                    ++references[reference_id].cov.bins_height[bin_number];
                }
            }
        }
    }
    matches_count = reads.size();

    float totalAb = 0.0;
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].reads_count > 0)
        {
            ++reference_count;
            matched_ref_length += references[i].length;
            references[i].abundance = float(references[i].reads_count * 100)/hits_count;
            totalAb += references[i].abundance/references[i].length;
        }
        else
        {
            references[i].abundance = 0.0;
        }
    }
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].reads_count > 0)
        {
            references[i].abundance = (references[i].abundance * 100) / (totalAb*references[i].length);
        }
    }

    totalAb = 0.0;
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].uniq_reads_count > 0)
        {
            references[i].uniq_abundance = float(references[i].uniq_reads_count * 100)/uniq_hits_count;
            totalAb += references[i].uniq_abundance/references[i].length;
        }
        else
        {
            references[i].uniq_abundance = 0.0;
        }
    }
    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].uniq_reads_count > 0)
        {
            references[i].uniq_abundance = (references[i].uniq_abundance * 100) / (totalAb*references[i].length);
        }

    }
}

//collect the sam files to process
inline void slimm::collect_bam_files()
{
    number_of_files = 1;
    if (options.is_directory)
    {
        _input_paths = get_bam_files_in_directory(options.input_path);
        number_of_files = length(_input_paths);
        if (options.verbose)
            std::cerr << number_of_files << " SAM/BAM Files found under the directory: " << options.input_path << "!\n";
    }
    else
    {
        if (is_file(toCString(options.input_path)))
            _input_paths.push_back(options.input_path);
        else
        {
            std::cerr << options.input_path << " is not a file use -d option for a directory.\n";
            exit(1);
        }
    }
}

float slimm::coverage_cut_off()
{
    if (_coverage_cut_off == 0.0 && options.cov_cut_off < 1.0)
    {
        std::vector<float> covs = {};
        covs.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                covs.push_back(references[i].cov_percent());
            }
        }
        _coverage_cut_off = get_quantile_cut_off<float>(covs, options.cov_cut_off);
    }
    return _coverage_cut_off;
}

float slimm::expected_coverage() const
{
    return float(avg_read_length * matches_count) / matched_ref_length;
}

inline void slimm::filter_alignments()
{
    uint32_t reference_count = length(references);
    for (uint32_t i=0; i < reference_count; ++i)
    {
        if (references[i].reads_count == 0)
            continue;
        if (references[i].cov_percent() >= coverage_cut_off() &&
            references[i].uniq_cov_percent() >= uniq_coverage_cut_off() && true)
        {
            valid_ref_ids.insert(i);
        }
        else
        {
            if(references[i].uniq_cov_percent() < uniq_coverage_cut_off())
            {
                ++failed_byUniqCov;
            }
            if(references[i].reads_count < options.min_reads)
            {
                ++failed_by_min_read;
            }
            if(references[i].cov_percent() < coverage_cut_off())
            {
                ++failed_byCov;
            }
        }
    }

    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        it->second.update(valid_ref_ids, references);
        if(it->second.is_uniq())
        {
            uint32_t reference_id = (it->second.targets[0]).reference_id;
            references[reference_id].uniq_reads_count2 += 1;
            uniq_matches_count2 += 1;
            uint32_t bin_number = (it->second.targets[0]).positions[0];
            ++references[reference_id].uniq_cov2.bins_height[bin_number];
        }
    }
}

// get taxonomic profiles from the sam/bam
inline void slimm::get_profiles()
{
    Timer<>  stop_watch;

    BamFileIn bam_file;
    BamHeader bam_header;

    std::cerr   << "\nReading " << current_file_index + 1 << " of " << number_of_files << " files ... ("
                << get_file_name(current_bam_file_path()) << ")\n"
                <<"=================================================================\n";

    if (read_bam_file(bam_file, bam_header, current_bam_file_path()))
    {
        //get average read length from a sample (size = 100K)
        avg_read_length = get_avg_read_length(bam_file, 100000);

        //if bin_width is not given use avg read length
        if (options.bin_width == 0)
            options.bin_width = avg_read_length;

        //reset the bam_file to the first recored by closing and reopening
        close(bam_file);

        read_bam_file(bam_file, bam_header, current_bam_file_path());

        StringSet<CharString>    contig_names = contigNames(context(bam_file));
        StringSet<uint32_t>      refLengths;
        refLengths = contigLengths(context(bam_file));

        uint32_t references_count = length(contig_names);
        references.resize(references_count);

        std::cerr<<"Intializing coverages for all reference genome ... ";
        // Intialize coverages for all genomes

        for (uint32_t i=0; i < references_count; ++i)
        {
            std::string accession = get_accession_id(contig_names[i]);
            uint32_t taxa_id = 0;
            auto ac_pos = db.ac__taxid.find(accession);
            if(ac_pos != db.ac__taxid.end())
            {
                taxa_id = ac_pos->second[0];
            }
            else
            {
                db.ac__taxid[accession] = std::vector<uint32_t>(LINAGE_LENGTH, 0);
            }
            reference_contig current_ref(accession, taxa_id, refLengths[i], options.bin_width);
            references[i] = current_ref;
        }
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

        std::cerr<<"Analysing alignments, reads and references ....... ";
        analyze_alignments(bam_file);
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
        if (hits_count == 0)
        {
            std::cerr << "[WARNING] No mapped reads found in BAM file!" << std::endl;
            return;
        }

        // Set the minimum reads to 10k-th of the total number of matched reads if not set by the user
        if (options.min_reads == 0)
          options.min_reads = 1 + ((matches_count - 1) / 10000);
        if (options.verbose)
            print_matches_stat();

        std::cerr   << "Filtering unlikely sequences ..................... ";
        filter_alignments();
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

        if (options.verbose)
            print_filter_stat();

        if (options.raw_output)
        {
            std::cerr<<"Writing features to a file ....................... ";
            write_raw_stat();
            std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
        }

        if (options.coverage_output)
        {
            std::cerr<<"Writing coverage profiles to a file ....................... ";
            write_coverage();
            std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;
        }

        std::cerr<<"Assigning reads to Least Common Ancestor (LCA) ... ";
        get_reads_lca_count();
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

        std::cerr<<"Writing taxnomic profile(s) ...................... ";
        write_abundance();
        if (options.verbose)
            std::cerr<<"\n.................................................. ";
        std::cerr<<"[" << stop_watch.lap() <<" secs]"  << std::endl;

        std::cerr<<"[Done!] File took " << stop_watch.elapsed() <<" secs to process.\n";
    }
}

inline void slimm::get_considered_ranks()
{
    if(options.rank == "all")
    {
        for(uint32_t i=8; i>0; --i)
            considered_ranks.push_back(static_cast<taxa_ranks>(i-1));
    }
    else if (options.rank == "superkingdom")
    {
        considered_ranks.push_back(to_taxa_ranks(options.rank));
    }
    else
    {
        considered_ranks.push_back(taxa_ranks(to_taxa_ranks(options.rank) + 1));
        considered_ranks.push_back(to_taxa_ranks(options.rank));
    }
}

inline uint32_t slimm::get_lca(std::set<uint32_t> const & ref_ids)
{
    uint32_t taxa_id = 1;
    for (uint32_t i=0; i<LINAGE_LENGTH; ++i)
    {
        std::set<uint32_t> level_taxa_set = {};
        for(auto ref_id : ref_ids)
        {
            taxa_id = db.ac__taxid[references[ref_id].accession][i];
            level_taxa_set.insert(taxa_id);
        }
        if(level_taxa_set.size() == 1)
            break;
    }
    return taxa_id;
}

inline void slimm::get_reads_lca_count()
{
    // put the non-unique read to upper taxa.
    for (auto it= reads.begin(); it != reads.end(); ++it)
    {
        size_t len = it->second.targets.size();
        if(len > 1)
        {
            uint32_t lca_taxa_id = 0;
            std::set<uint32_t> ref_ids = {};
            std::set<uint32_t> taxon_ids;
            size_t len = it->second.targets.size();
            for (size_t i=0; i < len; ++i)
            {
                uint32_t ref_id = (it->second.targets[i]).reference_id;
                ref_ids.insert(ref_id);
            }
            lca_taxa_id = get_lca(ref_ids);

            increment_or_initialize(taxon_id__read_count, lca_taxa_id, 1u);

            //add the contributing children references to the taxa
            taxon_id__children[lca_taxa_id].insert(ref_ids.begin(), ref_ids.end());
        }
    }

    //add the sum of read counts of children to all ancestors of the LCA // but get a copy first
    std::unordered_map <uint32_t, uint32_t> taxon_id__read_count_cp = taxon_id__read_count;
    uint32_t reciever_taxa_id = 0;
    for (auto t_id : taxon_id__read_count_cp)
    {
        // get the rank of the taxid
        taxa_ranks rnk = std::get<0>(db.taxid__name[t_id.first]);

        //get the first child and then the linage
        std::string first_child_acc = "";
        for (auto child : taxon_id__children.at(t_id.first))
        {
            first_child_acc = references[child].accession;
            break;
        }
        std::vector<uint32_t> linage = db.ac__taxid[first_child_acc];
        std::set<uint32_t> ref_ids = taxon_id__children[t_id.first];

        // add the read count to the uper ranks along the linage
        for (uint32_t j=rnk+1; j < LINAGE_LENGTH; ++j)
        {
            reciever_taxa_id = linage[j];
            increment_or_initialize(taxon_id__read_count, reciever_taxa_id, t_id.second);

            //add the contributing children references to the taxa
            taxon_id__children[reciever_taxa_id].insert(ref_ids.begin(), ref_ids.end());
        }
    }


    for (uint32_t i=0; i<length(references); ++i)
    {
        if (references[i].uniq_reads_count2 > 0)
        {
            std::vector<uint32_t> linage = db.ac__taxid[references[i].accession];
            std::set<uint32_t> ref_ids = taxon_id__children[linage[0]];
            for (uint32_t j=1; j<LINAGE_LENGTH; ++j)
            {
                reciever_taxa_id = linage[j];
                auto tid_pos = taxon_id__read_count.find(reciever_taxa_id);
                // If taxon_id already exists increment it
                if(tid_pos != taxon_id__read_count.end())
                    tid_pos->second += references[i].uniq_reads_count2;
                else
                    taxon_id__read_count[reciever_taxa_id] = references[i].uniq_reads_count2;

                //add the contributing children references to the taxa
                taxon_id__children[reciever_taxa_id].insert(i);
                taxon_id__children[reciever_taxa_id].insert(ref_ids.begin(), ref_ids.end());
            }
        }
    }
}

inline void slimm::print_filter_stat()
{
    std::cerr << "  " << length(valid_ref_ids) << " passed the threshould coverage.\n";
    std::cerr << "  " << failed_byCov << " ref's couldn't pass the coverage threshould.\n";
    std::cerr << "  " << failed_byUniqCov << " ref's couldn't pass the uniq coverage threshould.\n";
    std::cerr << "  uniquily matching reads increased from " << uniq_matches_count << " to " << uniq_matches_count2 <<"\n\n";
}

inline void slimm::print_matches_stat()
{
    std::cerr << "  "   << hits_count << " records processed." << std::endl;
    std::cerr << "    " << matches_count << " matching reads" << std::endl;
    std::cerr << "    " << uniq_matches_count << " uniquily matching reads"<< std::endl;
    std::cerr << "  references with reads = " << reference_count << std::endl;
    std::cerr << "  expected bins coverage = " << expected_coverage() <<std::endl;
    std::cerr << "  bins coverage cut-off = " << coverage_cut_off() << " (" << options.cov_cut_off <<" quantile)\n";
    std::cerr << "  uniq bins coverage cut-off = " << uniq_coverage_cut_off() << " (" << options.cov_cut_off <<" quantile)\n\n";
}

uint32_t slimm::min_reads()
{
    if (options.cov_cut_off == 1.0)
        _min_reads = 0;
    if (_min_reads == -1)
    {
        std::vector<int> counts = {};
        counts.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].reads_count > 0)
            {
                counts.push_back(references[i].reads_count);
            }
        }
        _min_reads = get_quantile_cut_off(counts, options.cov_cut_off);
    }
    return _min_reads;
}

uint32_t slimm::min_uniq_reads()
{
    if (options.cov_cut_off == 1.0)
        _min_uniq_reads = 0;
    if (_min_uniq_reads == -1)
    {
        std::vector<int> uniqCounts = {};
        uniqCounts.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                uniqCounts.push_back(references[i].uniq_reads_count);
            }
        }
        _min_uniq_reads = get_quantile_cut_off(uniqCounts, options.cov_cut_off);
    }
    return _min_uniq_reads;
}

float slimm::uniq_coverage_cut_off()
{
    if (_uniq_coverage_cut_off == 0.0 && options.cov_cut_off < 1.0)
    {
        std::vector<float> covs = {};
        covs.reserve(length(references));
        for (uint32_t i=0; i<length(references); ++i)
        {
            if (references[i].uniq_reads_count > 0)
            {
                covs.push_back(references[i].uniq_cov_percent());
            }
        }
        _uniq_coverage_cut_off = get_quantile_cut_off<float>(covs, options.cov_cut_off);
    }
    return _uniq_coverage_cut_off;
}

std::string slimm::get_lineage_string (taxa_ranks rank, std::vector<uint32_t> const & linage)
{
    std::string taxon_name = std::get<1>(db.taxid__name[linage[rank]]);
    if (taxon_name == "")
    {
        taxon_name =  "unknown_" + from_taxa_ranks(rank);
    }

    std::string linage_str = from_taxa_ranks_short(rank) + "__" + taxon_name;

    for (uint32_t i=rank+1; i < LINAGE_LENGTH; ++i)
    {
        taxon_name = std::get<1>(db.taxid__name[linage[i]]);
        if (taxon_name == "")
        {
            taxon_name =  "unknown_" + from_taxa_ranks(taxa_ranks(i));
        }
        linage_str = from_taxa_ranks_short(taxa_ranks(i)) + "__" + taxon_name + "|" + linage_str;
    }
    return linage_str;
}

std::string slimm::get_lineage_string (taxa_ranks rank, uint32_t const & taxa_id)
{
    std::vector<uint32_t> linage;
    if(taxa_id == 0)
    {
        linage.resize(LINAGE_LENGTH, 0);
    }
    else
    {
        std::string child_acc = "";
        for (auto child : taxon_id__children.at(taxa_id))
        {
            child_acc = references[child].accession;
            break;
        }
        linage = db.ac__taxid[child_acc];
    }
    return get_lineage_string(rank, linage);
}


inline void slimm::write_abundance()
{
    std::string abundunce_tsv_path = get_tsv_file_name(toCString(options.output_prefix), current_bam_file_path(), "_profile");
    std::ofstream abundunce_stream(abundunce_tsv_path);
    abundunce_stream << "taxa_level\ttaxa_id\tlinage\tabundance\tread_count\n";

    taxa_ranks rank = considered_ranks[1];
    taxa_ranks parent_rank = considered_ranks[0];

    // reserve the statics of un upper level
    std::unordered_map<uint32_t, float>     parent_abundance;
    std::unordered_map<uint32_t, uint32_t>  parent_reads_count;

    //get a hold of information at the upper taxon level
    for (auto t_id : taxon_id__read_count)
    {
        if (std::get<0>(db.taxid__name[t_id.first]) == parent_rank)
        {
            uint32_t genome_Length = 0;
            uint32_t children_count = 0;
            for (auto child : taxon_id__children.at(t_id.first))
            {
                genome_Length += references[child].length;
                ++children_count;
            }
            genome_Length = genome_Length/children_count;

            float abundance = float(t_id.second)/(matches_count) * 100;
            // New resolution into the unclassifieds
            parent_abundance[t_id.first] = abundance;
            parent_reads_count[t_id.first] = t_id.second;
        }
    }


    uint32_t    count = 0;
    uint32_t    faild_count = 0;
    uint32_t    sum_reads_count = 0.0;
    float       sum_abundunce = 0.0;

    std::unordered_map <uint32_t, float>    sum_abundunce_by_parent;
    std::unordered_map <uint32_t, uint32_t> sum_reads_count_by_parent;

    for (auto t_id : taxon_id__read_count)
    {
        if (std::get<0>(db.taxid__name[t_id.first]) == rank)
        {
            uint32_t genome_Length = 0;
            uint32_t children_count = 0;
            std::string child_acc = "";
            for (auto child : taxon_id__children.at(t_id.first))
            {
                genome_Length += references[child].length;
                child_acc = references[child].accession;
                ++children_count;
            }
            genome_Length = genome_Length/children_count;

            std::vector<uint32_t> linage = db.ac__taxid[child_acc];
            float cov = float(t_id.second * avg_read_length)/genome_Length;
            float abundance = float(t_id.second)/(matches_count) * 100;
            std::string candidate_name = std::get<1>(db.taxid__name[t_id.first]);

            // agregate the statstics of the children by parent
            uint32_t parent_tax_id = linage[parent_rank];
            increment_or_initialize (sum_abundunce_by_parent, parent_tax_id, abundance);
            increment_or_initialize (sum_reads_count_by_parent, parent_tax_id, t_id.second);
            if (abundance < options.abundance_cut_off || cov < coverage_cut_off() || candidate_name == "")
            {
                ++faild_count;
                continue;
            }
            std::string linage_str = get_lineage_string(rank, t_id.first);
            abundunce_stream << from_taxa_ranks(rank) << "\t" << t_id.first << "\t" << linage_str << "\t";
            abundunce_stream << abundance << "\t" << t_id.second << "\n";

            sum_abundunce += abundance;
            sum_reads_count += t_id.second;
            ++count;
        }
    }

    // unclassifieds with known parent
    for(auto ab_by_parent : sum_abundunce_by_parent)
    {
        uint32_t parent_taxid = ab_by_parent.first;
        float uncl_abundance = parent_abundance[parent_taxid] - sum_abundunce_by_parent[parent_taxid];
        uint32_t unc_read_count = parent_reads_count[parent_taxid] - sum_reads_count_by_parent[parent_taxid];
        std::string candidate_name = std::get<1>(db.taxid__name[parent_taxid]) + "_unclassified";
        if (uncl_abundance > options.abundance_cut_off && candidate_name != "_unclassified")
        {
            std::string linage_str = get_lineage_string(parent_rank, parent_taxid) + "|" + from_taxa_ranks_short(rank) + "__" + candidate_name;

            abundunce_stream << from_taxa_ranks(rank) << "\t" << parent_taxid << "*\t" << linage_str << "\t";
            abundunce_stream << uncl_abundance << "\t" << unc_read_count << "\n";
            sum_reads_count += unc_read_count;
            sum_abundunce += uncl_abundance;
        }
    }

    std::string linage_str = get_lineage_string(rank, 0);
    abundunce_stream << from_taxa_ranks(rank) << "\t" << "0*" << "\t" << linage_str << "\t";
    abundunce_stream << 100.0 - sum_abundunce << "\t" << matches_count - sum_reads_count << "\n";
    if (options.verbose)
    {
        std::cerr << "\n" << std::setw (4) << count << std::setw (15) << from_taxa_ranks(rank) <<" ("
        << faild_count <<" bellow cutoff i.e. "<< options.abundance_cut_off <<")";
    }

    abundunce_stream.close();
}


inline void slimm::write_coverage()
{
    std::string coverge_csv_path = get_tsv_file_name(options.output_prefix, current_bam_file_path(), "_coverge");
    std::string uniq_coverge_csv_path = get_tsv_file_name(options.output_prefix, current_bam_file_path(), "_uniq_coverge");
    std::string uniq_coverge2_csv_path = get_tsv_file_name(options.output_prefix, current_bam_file_path(), "_uniq_coverge2");

    std::ofstream coverge_stream(coverge_csv_path);
    std::ofstream uniq_coverge_stream(uniq_coverge_csv_path);
    std::ofstream uniq_coverge2_stream(uniq_coverge2_csv_path);


    for (auto valid_id : valid_ref_ids)
    {
        reference_contig current_ref = references[valid_id];
        coverge_stream  << current_ref.accession;
        uniq_coverge_stream  << current_ref.accession;
        uniq_coverge2_stream  << current_ref.accession;
        for (uint32_t b=0; b < current_ref.cov.number_of_bins; ++b)
        {
            coverge_stream  << "," << current_ref.cov.bins_height[b];
            uniq_coverge_stream  << "," << current_ref.uniq_cov.bins_height[b];
            uniq_coverge2_stream  << "," << current_ref.uniq_cov2.bins_height[b];
        }
        coverge_stream  << "\n" ;
        uniq_coverge_stream  << "\n";
        uniq_coverge2_stream  << "\n";
    }
    coverge_stream.close();
    uniq_coverge_stream.close();
    uniq_coverge2_stream.close();
}

inline void slimm::write_raw_stat()
{
    std::string raw_tsv_path = get_tsv_file_name(options.output_prefix, current_bam_file_path(), "_raw");
    std::ofstream features_stream(raw_tsv_path);

    features_stream <<"accesion\t"
                      "taxaid\t"
                      "name\t"
                      "reads_count\t"
                      "abundance\t"
                      "uniq1_abundance\t"
                      "uniq2_abundance\t"
                      "genome_length\t"
                      "uniq1_reads_count\t"
                      "uniq2_reads_count\t"

                      "bins_count\t"
                      "bins_count(>0)\t"
                      "uniq1_bins_count(>0)\t"
                      "uniq2_bins_count(>0)\t"

                      "coverage_depth\t"
                      "uniq1_coverage_depth\t"
                      "uniq2_coverage_depth\t"
                      "coverage(%)\t"
                      "uniq1_coverage(%)\t"
                      "uniq2_coverage(%)\n";

    for (uint32_t i=0; i < length(references); ++i)
    {
        reference_contig current_ref = references[i];
        std::string candidate_name = std::get<1>(db.taxid__name[current_ref.taxa_id]);
        if (candidate_name == "")
            candidate_name = "no_name_found";
        features_stream   << current_ref.accession << "\t"
                          << current_ref.taxa_id << "\t"
                          << candidate_name << "\t"
                          << current_ref.reads_count << "\t"
                          << current_ref.abundance << "\t"
                          << current_ref.uniq_abundance << "\t"
                          << current_ref.uniq_abundance2 << "\t"
                          << current_ref.length << "\t"
                          << current_ref.uniq_reads_count << "\t"
                          << current_ref.uniq_reads_count2 << "\t"

                          << current_ref.cov.number_of_bins << "\t"
                          << current_ref.cov.none_zero_bin_count() << "\t"
                          << current_ref.uniq_cov.none_zero_bin_count() << "\t"
                          << current_ref.uniq_cov2.none_zero_bin_count() << "\t"

                          << current_ref.cov_depth() << "\t"

                          << current_ref.uniq_cov_depth() << "\t"
                          << current_ref.uniq_cov_depth2() << "\t"

                          << current_ref.cov_percent() << "\t"
                          << current_ref.uniq_cov_percent() << "\t"
                          << current_ref.uniq_cov_percent2() << "\n";
    }
    features_stream.close();
}


inline int get_taxonomic_profile(arg_options & options)
{
    // slimm object
    Timer<>  stop_watch;
    uint32_t total_hits_count = 0;
    slimm slimm1(options);
    for (uint32_t n=0; n < slimm1.number_of_files; ++n)
    {
        slimm1.reset();
        slimm1.current_file_index = n;
        slimm1.get_profiles();
        total_hits_count += slimm1.hits_count;
    }

    std::string output_directory = get_directory(options.output_prefix);

    std::cerr << "\n*****************************************************************\n";
    std::cerr << total_hits_count << " SAM/BAM alignment records are proccessed.\n";
    std::cerr << "Taxonomic profiles are written to: \n   " << output_directory <<"\n";
    std::cerr << "Total time elapsed: " << stop_watch.elapsed() <<" secs\n";

    return 0;
}


#endif /* SLIMM_H */
