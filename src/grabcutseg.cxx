#include <common.h>
#include <gmm.h>
#include <graphcuts.h>
#include <utility.h>
#include <loadfiles.h>
#include <query.h>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
     std::string input_image_file, init_image_file, alpha_file, tlink_image_file, nlink_image_file, mask_file, gmm_fg_file, gmm_bg_file, priorfg_file, priorbg_file, fg_la