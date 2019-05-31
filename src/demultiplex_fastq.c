/*
 * Simple, fast and memory efficient demultiplexer for FASTQ sequencing files
 * written by Andy Hauser <andreas.hauser@LMU.de>
 */

#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <zlib.h>

#include <execinfo.h>
#include <signal.h>

#include "kseq.h"


KSEQ_INIT(gzFile, gzread)
void demultiplex_debug_catch_signal(int sig_num)
{
    fprintf (stderr, "\n\n-------------------------8<-----------------------\nExiting on error!\n");
    fprintf (stderr, "Signal %d received\n",sig_num);
    perror("ERROR (can be bogus)");
    fprintf(stderr, "Backtrace:");
    void *buffer[30];
    int nptrs = backtrace(buffer, 30);
    backtrace_symbols_fd(buffer, nptrs, 2);
    fprintf (stderr, "------------------------->8-----------------------\n\n"
        "Send the binary program that caused this error and the coredump (ls core.*).\n"
        "Or send the backtrace:"
        "\n$ gdb -ex=bt --batch PROGRAMM_NAME CORE_FILE\n"
        "If there is no core file, enable coredumps in your shell and run again:\n"
        "$ ulimit -c unlimited\n\n");
    fprintf (stderr, "Please report this to Andy Hauser <Andreas.Hauser@lmu.de>.\n");

  //exit(sig_num);
  abort();
}



/* See: https://stackoverflow.com/questions/9210528/split-string-with-delimiters-in-c */
char** str_split(char* a_str, const char a_delim, size_t *count_res)
{
  char** result    = 0;
  size_t count     = 0;
  char* tmp        = a_str;
  char* last_comma = 0;
  char delim[2];
  delim[0] = a_delim;
  delim[1] = 0;

  while (*tmp)
  {
    if (a_delim == *tmp)
    {
      count++;
      last_comma = tmp;
    }
    tmp++;
  }

  count += last_comma < (a_str + strlen(a_str) - 1);

  *count_res = count;

  count++;

  result = calloc(count + 1, sizeof(char*));

  if (result)
  {
    size_t idx  = 0;
    char* token = strtok(a_str, delim);

    while (token)
    {
      assert(idx < count);
      *(result + idx++) = strndup(token, FILENAME_MAX);
      token = strtok(0, delim);
    }
    assert(idx == count - 1);
    *(result + idx) = 0;
  }

  return result;
}

char *prefix(char* path)
{
  char* name = basename(strndup(path, FILENAME_MAX));
  char* suffix = strrchr(name, '.');
  if(suffix)
  {
    if(strcmp(suffix, "gz"))
    {
      *suffix = '\0';
      return prefix(name); // XXX free name
    }
    *suffix = '\0';
  }
  return name;
}

void print_indexed_seq(FILE *fout1, kseq_t *seq1, kseq_t *seqi, kseq_t *seqj)
{
  if (seq1->comment.l)
    if(seqj)
      fprintf(fout1, "@%s %s %s:%s\n", seq1->name.s, seq1->comment.s, seqi->seq.s, seqj->seq.s);
    else
      fprintf(fout1, "@%s %s %s\n", seq1->name.s, seq1->comment.s, seqi->seq.s);
  else
    if(seqj)
      fprintf(fout1, "@%s %s:%s\n", seq1->name.s, seqi->seq.s, seqj->seq.s);
    else
      fprintf(fout1, "@%s %s\n", seq1->name.s, seqi->seq.s);

  fprintf(fout1, "%s\n", seq1->seq.s);

  if (seq1->qual.l)
    fprintf(fout1, "+\n%s\n", seq1->qual.s);
  else
    fprintf(fout1, "+\n\n");
}

/* Flag set by ‘--verbose’. */
int verbose_flag;

int main (int argc, char **argv)
{
  char *fastq1 = NULL, *fastq2 = NULL, *fastqi = NULL, *fastqj = NULL, *prefix_out = NULL, *prefix_out2 = NULL;
  char **barcodes = NULL;
  char **barcodesj = NULL;
  FILE **barcode_files1 = NULL;
  FILE **barcode_files2 = NULL;
  size_t n_barcodes = 0;
  gzFile fq1 = NULL, fq2 = NULL, fqi = NULL, fqj = NULL;
  kseq_t *seq1, *seq2 = NULL, *seqi, *seqj = NULL;
  int l1,l2,li,lj;
  /* Flag set by ‘--no_other’. */
  int no_other = 0;


  int c;

  struct sigaction handler;
  handler.sa_handler = demultiplex_debug_catch_signal;
  sigemptyset(&handler.sa_mask);
  handler.sa_flags = 0;

  sigaction(SIGFPE, &handler, NULL);
  sigaction(SIGSEGV, &handler, NULL);
  sigaction(SIGBUS, &handler, NULL);


  while (1)
  {
    static struct option long_options[] =
    {
      {"verbose",    no_argument,       &verbose_flag, 1},
      {"no_other",   no_argument,       0, 'n'},
      {"r1",         required_argument, 0, '1'},
      {"r2",         required_argument, 0, '2'},
      {"i1",         required_argument, 0, 'i'},
      {"i2",         required_argument, 0, 'j'},
      {"barcodes",   required_argument, 0, 'b'},
      {"prefix",     required_argument, 0, 'p'},
      {"prefix2",    required_argument, 0, 'q'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "v1:2:i:j:b:p:q:n", long_options, &option_index);

    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;
        //printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
        break;

      case 'i':
        //printf ("option -i with value `%s'\n", optarg);
        fastqi = optarg;
        break;

      case 'j':
        //printf ("option -i with value `%s'\n", optarg);
        fastqj = optarg;
        break;

      case '1':
        //printf ("option -1 with value `%s'\n", optarg);
        fastq1 = optarg;
        break;

      case '2':
        //printf ("option -2 with value `%s'\n", optarg);
        fastq2 = optarg;
        break;

      case 'n':
        no_other = 1;
        break;

      case 'p':
        //printf ("option -p with value `%s'\n", optarg);
        prefix_out = optarg;
        break;

      case 'q':
        //printf ("option -p with value `%s'\n", optarg);
        prefix_out2 = optarg;
        break;

      case 'b':
        //printf ("option -b with value `%s'\n", optarg);
        barcodes = str_split(optarg,',', &n_barcodes);
        barcodesj = calloc(n_barcodes + 1, sizeof(char*));
        for(size_t i = 0; barcodes[i] != NULL; i++)
        {
          char* barcode = barcodes[i];
          char* second = strchr(barcode, ':');
          if(second)
          {
            *second = '\0';
            barcodesj[i] = second + 1;
          }
        }

        break;

      case '?':
        break;

      default:
        fprintf(stderr, "unknown argument: %c = %s\n", c, optarg);
        return EXIT_FAILURE;
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    printf ("non-option ARGV-elements: ");
    while (optind < argc)
      printf ("%s ", argv[optind++]);
    putchar ('\n');
    return EXIT_FAILURE;
  }

  if (fastq1 == NULL || barcodes == NULL)
  {
    fprintf(stderr, "Usage: %s [-n] [-p PREFIX] --r1 FASTQ1 [--r2 FASTQ2] --i1 FASTQ_INDEX [--i2 FASTQ_INDEX2] -b BARCODE1,BARCODE2[,...]\n", argv[0]);
    fprintf(stderr, "  --no_other,-n	Do not output non matching barcodes into OTHER file\n"
                    "  --prefix,-p PREFIX	prefix output filenames with PREFIX\n"
                    "  --r1 FASTQ1 	fastq file with reads\n"
                    "  --r2 FASTQ2	paired fastq file with reads\n"
                    "  --i1 FASTQ_INDEX	fist index read\n"
                    "  --i2 FASTQ_INDEX2	second index read\n"
                    "  --barcodes, -b BARCODE1,BARCODE2,...	barcodes, comma separated\n"
                    "\n	written by Andreas.Hauser@LMU.de\n");
    return EXIT_FAILURE;
  }

  prefix_out = prefix(fastq1);
  if(fastq2)
    prefix_out2 = prefix(fastq2);

  size_t *barcode_len = calloc(n_barcodes + 1, sizeof(size_t));
  size_t *barcodej_len = calloc(n_barcodes + 1, sizeof(size_t));
  barcode_files1 = calloc(n_barcodes + 1, sizeof(FILE *));
  barcode_files2 = calloc(n_barcodes + 1, sizeof(FILE *));
  char *filename1 = calloc(FILENAME_MAX,sizeof(char));
  char *filename2 = calloc(FILENAME_MAX,sizeof(char));
  int i;
  if(!no_other)
  {
    barcodes[n_barcodes] = "OTHER"; // pseudobarcode, at last entry
    barcodesj[n_barcodes] = "OTHER"; // pseudobarcode, at last entry
  }
  for(i = 0; barcodes[i] != NULL; i++)
  {
    barcode_len[i]  = strnlen(barcodes[i],FILENAME_MAX);
    strncpy(filename1,prefix_out,FILENAME_MAX);
    strncat(filename1,"_",FILENAME_MAX);
    strncat(filename1,barcodes[i],FILENAME_MAX);

    if(barcodesj[i])
    {
      barcodej_len[i] = strnlen(barcodesj[i],FILENAME_MAX);
      strncat(filename1,":",FILENAME_MAX);
      strncat(filename1,barcodesj[i],FILENAME_MAX);
    }
    strncat(filename1,".fastq",FILENAME_MAX);
    barcode_files1[i]  = fopen(filename1, "a");
    if(!barcode_files1[i])
    {
      perror("Opening file");
      exit(1);
    }

    if(fastq2)
    {
      strncpy(filename2,prefix_out2,FILENAME_MAX);
      strncat(filename2,"_",FILENAME_MAX);
      strncat(filename2,barcodes[i],FILENAME_MAX);
      if(barcodesj[i])
      {
        strncat(filename2,":",FILENAME_MAX);
        strncat(filename2,barcodesj[i],FILENAME_MAX);
      }
      strncat(filename2,".fastq",FILENAME_MAX);
      barcode_files2[i] = fopen(filename2, "a");
    }
  }

  fq1 = gzopen(fastq1, "r"); if(!fq1) perror(fastq1); seq1 = kseq_init(fq1); if(!seq1) perror(fastq1);
  fqi = gzopen(fastqi, "r"); if(!fqi) perror(fastqi); seqi = kseq_init(fqi); if(!seqi) perror(fastqi);

  if(fastq2)
  {
    fq2 = gzopen(fastq2, "r"); seq2 = kseq_init(fq2);
  }

  if(fastqj)
  {
    fqj = gzopen(fastqj, "r");
    seqj = kseq_init(fqj);
  }

  /*
  for(int i = 0; barcodes[i] != NULL; i++)
    if(barcodesj[i])
      fprintf(stderr, "%s %s\n", barcodes[i], barcodesj[i]);
    else
      fprintf(stderr, "%s\n", barcodes[i]);
   */

  while ( (l1 = kseq_read(seq1)) >= 0 &&
          (li = kseq_read(seqi)) >= 0 )
  {
    if(fastq2 && !((l2 = kseq_read(seq2)) >= 0))
      break;
    if(fastqj && !((lj = kseq_read(seqj)) >= 0))
      break;
    FILE *fout1, *fout2;
    if((strncmp(seq1->name.s, seqi->name.s,200) != 0) ||
       (fq2 && strncmp(seq1->name.s, seq2->name.s,200) != 0)
       )
    {
      fprintf(stderr, "name mismatch:\n");
      fprintf(stderr, "seq 1 name: %s\n", seq1->name.s);
      fprintf(stderr, "seq 2 name: %s\n", seq2->name.s);
      fprintf(stderr, "index 1 name: %s\n", seqi->name.s);
    }

    int i = 0;
    if(fastqj)
      for(i = 0; barcodes[i] != NULL; i++)
      {
        if(strncmp(seqi->seq.s, barcodes[i],  barcode_len[i])  == 0 &&
           strncmp(seqj->seq.s, barcodesj[i], barcodej_len[i]) == 0 )
          break;
      }
    else
      for(i = 0; barcodes[i] != NULL; i++)
        if(strncmp(seqi->seq.s, barcodes[i], barcode_len[i]) == 0)
          break;

    if(barcodes[i] == NULL)
    {
      if(no_other == 0)
        i--; // use pseudobarcode at last entry in barcode mapping
      else
        continue;
    }

    fout1   = barcode_files1[i];
    print_indexed_seq(fout1, seq1, seqi, seqj);
    if(fastq2)
    {
      fout2    = barcode_files2[i];
      print_indexed_seq(fout2, seq2, seqi, seqj);
    }
  }

  kseq_destroy(seq1);
  kseq_destroy(seqi);
  gzclose(fq1);
  gzclose(fqi);
  if(fq2)
  {
    kseq_destroy(seq2);
    gzclose(fq1);
  }
  return 0;
}
