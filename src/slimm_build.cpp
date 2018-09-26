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
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cstdlib>


#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include "misc.hpp"
#include "file_helper.hpp"

using namespace seqan;

// ----------------------------------------------------------------------------
// Class arg_options
// ----------------------------------------------------------------------------
struct arg_options
{
    std::string                  fasta_path;
    std::string                  nodes_path;
    std::string                  names_path;
    std::string                  output_path;

    arg_options() : fasta_path(),
                    nodes_path(),
                    names_path(),
                    output_path("slimm_db.sldb"){}
};

// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
void setupArgumentParser(ArgumentParser & parser, arg_options const & options)
{
    // Setup ArgumentParser.
    setAppName(parser, "slimm_build");
    setShortDescription(parser, "gets a reduced taxonomic information given a multi-fasta file using accession numbers");
    setCategory(parser, "Metagenomics");

    setDateAndVersion(parser);
    setDescription(parser);
    // Define usage line and long description.
    addUsageLine(parser, "-nm \"\\fINAMES.dmp\\fP\" -nd \"\\fINODES.dmp\\fP\" -o \"\\fISLIMM.sldb\\fP\" [\\fIOPTIONS\\fP] \"\\fIFASTA\\fP]");

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "FASTA FILE"));
    setValidValues(parser, 0, SeqFileIn::getFileExtensions());
    setHelpText(parser, 0, "A multi-fasta file used as a reference for mapping");

    // The output file argument.
    addOption(parser, ArgParseOption("o", "output-file", "The path to the output file (default slimm_db.sldb)",
                                     ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "output-file", ".sldb");
    setDefaultValue(parser, "output-file", options.output_path);

    addOption(parser, ArgParseOption("nm", "names", "NCBI's names.dmp file which contains the mapping of taxaid to name",
                             ArgParseArgument::INPUT_FILE));
    setRequired(parser, "names");

    addOption(parser, ArgParseOption("nd", "nodes", "NCBI's nodes.dmp file which contains the taxonomic tree.",
                             ArgParseArgument::INPUT_FILE));
    setRequired(parser, "nodes");
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------
ArgumentParser::ParseResult
parseCommandLine(ArgumentParser & parser, arg_options & options, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    getArgumentValue(options.fasta_path, parser, 0);

    getOptionValue(options.names_path, parser, "names");
    getOptionValue(options.nodes_path, parser, "nodes");

    if (isSet(parser, "output-file"))
        getOptionValue(options.output_path, parser, "output-file");

    return ArgumentParser::PARSE_OK;
}


// --------------------------------------------------------------------------
// Function get_accession_numbers()
// --------------------------------------------------------------------------
inline void get_accession_numbers(std::vector<std::string> & accessions, arg_options const & options)
{

    std::set<std::string> accession_set;
    std::cerr <<"[MSG] getting accessions numbers from fasta file ...\n";
    CharString id;
    IupacString seq;

    SeqFileIn fasta_file;
    if (!open(fasta_file, toCString(options.fasta_path)))
    {
        CharString msg = "Unable to open contigs File: ";
        append (msg, options.fasta_path);
        throw toCString(msg);
    }
    while(!atEnd(fasta_file))
    {

        readRecord(id, seq, fasta_file);
        accession_set.insert(get_accession_id(id));
    }
    close(fasta_file);

    // fill the std::vector from the set
    for(auto ac_it=accession_set.begin(); ac_it != accession_set.end(); ++ac_it)
    {
        accessions.push_back(*ac_it);
    }
}


// --------------------------------------------------------------------------
// Function print_missed_accessions()
// --------------------------------------------------------------------------
inline void print_missed_accessions(std::vector<std::string> & accessions,
                                    arg_options const & options)
{
    uint32_t count = 3;
    uint32_t dot_pos = options.output_path.size() - 4;
    std::string missed_acc_path = options.output_path.substr(0, dot_pos) + "missed";
    std::ofstream missed_acc_stream(missed_acc_path);

    std::cerr <<"[WARNING!] "<< accessions.size() <<" accessions (";
    for(auto ac_it=accessions.begin(); count > 0 && ac_it != accessions.end(); --count, ++ac_it)
        std::cerr << *ac_it << ", ";
    std::cerr <<"...) were not mapped to taxaid.\n";

    for(auto ac_it=accessions.begin(); ac_it != accessions.end(); ++ac_it)
        missed_acc_stream << *ac_it << "\n";
    missed_acc_stream.close();

    std::cerr <<"[WARNING!] Take a look at "<< missed_acc_path << " file for a complete list.\n";
}

// --------------------------------------------------------------------------
// Function get_taxid_from_xml_string()
// --------------------------------------------------------------------------
inline int get_taxid_from_xml_string(std::string & xml_str)
{
    int tax_id = -1;
    std::string const taxid_identifier = "<Item Name=\"TaxId\" Type=\"Integer\">";
    size_t tid_index = xml_str.find(taxid_identifier, 0);
    if (tid_index < xml_str.length())
    {
        size_t tid_begin = tid_index + taxid_identifier.length();
        size_t tid_end = xml_str.find('<', tid_begin);
        std::string tid_str = xml_str.substr(tid_begin, tid_end-tid_begin);
        tax_id = atoi( tid_str.c_str() ) ;
    }
    return tax_id;
}
// --------------------------------------------------------------------------
// Function download_xml_files()
// --------------------------------------------------------------------------
inline void download_xml_files(std::vector<std::string> & accessions, uint32_t const & start, uint32_t const & end)
{
    uint32_t const batch_size = end - start;
    uint32_t const num_threads = 20;
    uint32_t const max_download_attempts = 20;
    uint32_t chunk_size = batch_size/num_threads;
    if (chunk_size * num_threads != batch_size)
        chunk_size++;

    std::string const url_base = "\"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=";

    std::vector<std::future<void>> tasks;

    for (uint32_t task_no = 0; task_no < num_threads; ++task_no)
    {
        tasks.emplace_back(std::async([=, &accessions, &url_base] {
            for (uint32_t i = task_no*chunk_size; i < batch_size && i < (task_no +1)*chunk_size; ++i)
            {
                std::string cmd = "curl -s " + url_base + accessions[i+start] +"\" > " + accessions[i+start] +".xml";
                int download_attempts = 0;
                int download_res = 1;
                while (download_attempts < max_download_attempts && download_res == 1)
                {
                    download_res = std::system(cmd.c_str());
                }
            }
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }
}

inline void remove_xml_files(std::vector<std::string> & accessions, uint32_t const & start, uint32_t const & end)
{
    uint32_t const batch_size = end - start;
    uint32_t const num_threads = 20;
    uint32_t chunk_size = batch_size/num_threads;
    if (chunk_size * num_threads != batch_size)
        chunk_size++;

    std::vector<std::future<void>> tasks;

    for (uint32_t task_no = 0; task_no < num_threads; ++task_no)
    {
        tasks.emplace_back(std::async([=, &accessions] {
            for (uint32_t i = task_no*chunk_size; i < batch_size && i < (task_no +1)*chunk_size; ++i)
            {
                std::string cmd = "rm -rf " + accessions[i+start] +".xml";
                std::system(cmd.c_str());
            }
        }));
    }
    for (auto &&task : tasks)
    {
        task.get();
    }
}


// --------------------------------------------------------------------------
// Function get_taxid_from_accession()
// --------------------------------------------------------------------------
inline void get_taxid_from_accession(slimm_database & slimm_db,
                                     std::vector<std::string> & accessions,
                                     std::vector<std::string> & missed_accessions,
                                     arg_options const & options,
                                     uint32_t const & start,
                                     uint32_t const & end)
{
    for(uint32_t i=start; i<end; ++i)
    {
        std::string curr_acc = accessions[i];
        std::ifstream info_xml_stream(curr_acc + ".xml");
        if (info_xml_stream)
        {
            std::string xml_str((std::istreambuf_iterator<char>(info_xml_stream)),
                            std::istreambuf_iterator<char>());
            int taxa_id = get_taxid_from_xml_string(xml_str);
            if (taxa_id != -1)
            {
                slimm_db.ac__taxid[curr_acc] = std::vector<uint32_t>(LINAGE_LENGTH, 0);
                slimm_db.ac__taxid[curr_acc][0] = taxa_id;
            }
            else
            {
                std::cerr << "[" << curr_acc << "] BAD xml file doenloaded from entrez/eutils!" << std::endl;
                missed_accessions.push_back(curr_acc);
            }
        }
        else
        {
            std::cerr << "[ERR] " << curr_acc << ".xml not found. unable to use the entrez/eutils online service!" << std::endl;
            missed_accessions.push_back(curr_acc);
        }
        info_xml_stream.close();
    }
}

// --------------------------------------------------------------------------
// Function fill_name_taxid_linage()
// --------------------------------------------------------------------------
inline void fill_name_taxid_linage(slimm_database & slimm_db, arg_options const & options)
{
    std::unordered_map<uint32_t, std::tuple<taxa_ranks, uint32_t> > taxid__parent;
    std::unordered_map<uint32_t, std::string>                       taxid__name;

    std::ifstream taxid__parent_stream(options.nodes_path);
    std::ifstream taxid__name_stream(options.names_path);

    uint32_t taxid=0, parent_taxid = 0;
    std::string line, rank, ignore, name;

    while(std::getline(taxid__parent_stream, line))
    {
        std::stringstream   linestream(line);
        linestream >> taxid;// first column is taxid
        std::getline(linestream, ignore, '\t'); //skip |
        std::getline(linestream, ignore, '\t'); //skip |
        linestream >> parent_taxid; // third column is parent_taxid
        std::getline(linestream, ignore, '\t'); //skip |
        std::getline(linestream, ignore, '\t'); //skip |
        std::getline(linestream, rank, '\t'); //fifth column is rank
        taxid__parent[taxid] = std::make_tuple(to_taxa_ranks(rank), parent_taxid);
    }

    taxid__parent_stream.close();

    while(std::getline(taxid__name_stream, line))
    {
        auto pos = line.find("scientific name", 0);
        if(pos != std::string::npos)
        {
            std::stringstream   linestream(line);
            linestream >> taxid;// first column is taxid
            std::getline(linestream, ignore, '\t'); //skip |
            std::getline(linestream, ignore, '\t'); //skip |
            std::getline(linestream, name, '\t'); // 2nd column is name
            taxid__name[taxid] = name;
        }
    }
    taxid__name_stream.close();

    std::cerr <<"[MSG] getting taxonomic linages and resolving names ...\n";
    for(auto ac__taxid_it=slimm_db.ac__taxid.begin(); ac__taxid_it != slimm_db.ac__taxid.end(); ++ac__taxid_it)
    {
        uint32_t tid = ac__taxid_it->second[0];
        slimm_db.taxid__name[tid] = std::make_tuple(strain_lv, taxid__name[tid]);

        while (tid != 1)
        {
            auto tid_pos = taxid__parent.find(tid);
            if(tid_pos == taxid__parent.end())
                break;

            taxa_ranks current_rank = std::get<0>(tid_pos->second);
            if (current_rank >= species_lv && current_rank <= superkingdom_lv)
            {
                ac__taxid_it->second[current_rank] = tid;
                slimm_db.taxid__name[tid] = std::make_tuple(current_rank, taxid__name[tid]);
            }
            tid = std::get<1>(tid_pos->second);
        }
    }
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.
int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    arg_options options;
    setupArgumentParser(parser, options);

    ArgumentParser::ParseResult res = parseCommandLine(parser, options, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // get the accession numbers from the fasta file
    std::vector<std::string> accessions;
    std::vector<std::string> missed_accessions;
    get_accession_numbers(accessions, options);
    slimm_database slimm_db;

    uint32_t const batch_download = 200;
    uint32_t num_partitions = accessions.size()/batch_download;
    if (num_partitions * batch_download != accessions.size())
        num_partitions ++;

    std::cout << "[MSG] " << accessions.size()  << " accessions found! -- organized into " << num_partitions << " partitions!" << std::endl;

    for (uint32_t batch_no = 0; batch_no < num_partitions; ++batch_no)
    {

        uint32_t start = batch_no * batch_download;
        uint32_t end = ((batch_no+1) * batch_download > accessions.size())? accessions.size() : (batch_no+1) * batch_download;

        std::cerr <<"[MSG] Downloading accession xmls  \t" << batch_no+1 << "/" << num_partitions <<"\t..." << std::endl;
        download_xml_files(accessions, start, end);
        std::cerr <<"[MSG] mapping accessions to taxon id \t..." << std::endl;
        get_taxid_from_accession(slimm_db, accessions, missed_accessions, options, start, end);

        remove_xml_files(accessions, start, end);
    }

    // some accessions are still not mapped
    if(missed_accessions.size() > 0)
        print_missed_accessions(missed_accessions, options);

    std::cerr <<"[MSG] loading nodes and names mappings from files ..." << std::endl;
    fill_name_taxid_linage(slimm_db, options);
    std::cerr <<"[MSG] saving the SLIMM database to file ..." << std::endl;
    save_slimm_database(slimm_db, options.output_path);

//    std::vector<uint32_t> tids = slimm_db.ac__taxid["NZ_CP009257.1"];
//    std::cout << "ACC: NC_004578.1 \t taxa id: " << tids[0] << "\n";
//    for (uint32_t i=0; i<tids.size(); ++i)
//    {
//        std::string r = from_taxa_ranks(static_cast<taxa_ranks>(i));
//        std::cout << r << "\t" << from_taxa_ranks(std::get<0>(slimm_db.taxid__name[tids[i]])) << "\t" << tids[i] << "\t" << std::get<1>(slimm_db.taxid__name[tids[i]])  << "\n";
//    }

    return 0;
}
