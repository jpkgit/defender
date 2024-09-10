#include "baseline_alerts.h"
#include "utilities.h"

volatile bool do_exit = false;

FILE *outfile = NULL;
volatile uint32_t byte_count = 0;
volatile unsigned int sweep_count = 0;

struct timeval time_start;
struct timeval t_start;

static hackrf_device *device = NULL;
bool amp = false;
uint32_t amp_enable;

bool antenna = false;
uint32_t antenna_enable;

bool timestamp_normalized = false;
bool binary_output = false;
bool one_shot = false;
bool finite_mode = false;
volatile bool sweep_started = false;

int fftSize = 20;
double fft_bin_width;
fftwf_complex *fftwIn = NULL;
fftwf_complex *fftwOut = NULL;
fftwf_plan fftwPlan = NULL;
fftwf_complex *ifftwIn = NULL;
fftwf_complex *ifftwOut = NULL;
fftwf_plan ifftwPlan = NULL;
uint32_t ifft_idx = 0;
float *pwr;
float *window;
struct timeval usb_transfer_time;

/* New stuff*/
float baseline[BASELINE_SIZE];
float saved_baseline[BASELINE_SIZE];
bool baseline_saved = false;
bool b_quit = false;
pthread_mutex_t mutex;


static float TimevalDiff(const struct timeval *a, const struct timeval *b)
{
	return (a->tv_sec - b->tv_sec) + 1e-6f * (a->tv_usec - b->tv_usec);
}

int save_baseline()
{		
	fprintf(stderr, "Saving %u size baseline to file %s\n", BASELINE_SIZE, "/tmp/baseline.txt");
	
	int index = 0;
	for (index = 0; index < BASELINE_SIZE;index++)
	{
		saved_baseline[index] = -80.0;
 	}
	
	pthread_mutex_lock(&mutex);

	int j = 0;	
    for (j = 0; j < BASELINE_SIZE; j++) 
	{
        saved_baseline[j] = baseline[j];
    }
	
	baseline_saved = true;	
	pthread_mutex_unlock(&mutex);
    FILE *file = fopen("/tmp/baseline.txt", "w");

    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Write the array to the file
	int i = 0;

    for (i = 0; i < BASELINE_SIZE; i++) 
	{
        fprintf(file, "%f\n", saved_baseline[i]);
    }

    // Close the file
    fclose(file);
	fprintf(stderr, "Saved %u size baseline to file %s\n", BASELINE_SIZE, "/tmp/baseline.txt");
	return 0;	
}

int load_baseline()
 {
    FILE *file = fopen("/tmp/baseline.txt", "r");

    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }

	fprintf(stderr, "Loading %u size baseline from file %s\n", BASELINE_SIZE, "/tmp/baseline.txt");
	pthread_mutex_lock(&mutex);

	int i = 0;
    for (i = 0; i < BASELINE_SIZE; i++) 
	{
        saved_baseline[i] = -80.0;
    }
    
    char line[MAX_LINE_LENGTH];
	int index = 0;

    while (fgets(line, sizeof(line), file) != NULL) 
	{        
        // Convert line to float and store in array
        saved_baseline[index] = strtof(line, NULL);
        index++;
    }

	pthread_mutex_unlock(&mutex);
	fprintf(stderr, "Loaded %d size baseline from file %s\n", index, "/tmp/baseline.txt");
    fclose(file);
    return index;
}

int process_command(int command_key)
{
	switch (command_key)
	{
	case 'q':
		do_exit = true;
		break;
	case 's':
		save_baseline();
		break;
	case 'l':
		load_baseline();
		break;		
	case 'i':
		threshold++;
		fprintf(stderr, "Increased threshold to: [%d]\n", threshold);
		break;		
	case 'd':
		threshold--;
		fprintf(stderr, "Decreased threshold to: [%d]\n", threshold);
		break;		
	case 'a':
		average++;
		fprintf(stderr, "Increased average to: [%d]\n", average);			
		break;		
	case 'z':
		average--;
		fprintf(stderr, "Decreased average to: [%d]\n", average);	
		break;						
	case 'o':
		baseline_alert_offset++;
		fprintf(stderr, "Increased baseline alert offset to: [%d]\n", baseline_alert_offset);	
		break;		
	case 'p':
		baseline_alert_offset--;
		fprintf(stderr, "Decreased baseline alert offset to: [%d]\n", baseline_alert_offset);	
		break;	

	default:
		break;
	}

	return 0;
}

float logPower(fftwf_complex in, float scale)
{
	float re = in[0] * scale;
	float im = in[1] * scale;
	float magsq = re * re + im * im;
	return (float)(log2(magsq) * 10.0f / log2(10.0f));
}

int rx_callback(hackrf_transfer *transfer)
{
	int8_t *buf;
	uint8_t *ubuf;
	unsigned int frequency; /* in Hz */	
	int i, j;
	struct tm *fft_time;
	char time_str[50];	

	if (NULL == outfile)
	{
		return -1;
	}

	if (do_exit)
	{
		return 0;
	}

	// happens only once with timestamp_normalized == true
	if ((usb_transfer_time.tv_sec == 0 && usb_transfer_time.tv_usec == 0) ||
		timestamp_normalized == false)
	{
		// set the timestamp for the first sweep
		gettimeofday(&usb_transfer_time, NULL);
	}

	byte_count += transfer->valid_length;
	buf = (int8_t *)transfer->buffer;	

	for (j = 0; j < BLOCKS_PER_TRANSFER; j++)
	{
		ubuf = (uint8_t *)buf;

		if (ubuf[0] == 0x7F && ubuf[1] == 0x7F)
		{
			frequency = ((uint64_t)(ubuf[9]) << 56) |
						((uint64_t)(ubuf[8]) << 48) |
						((uint64_t)(ubuf[7]) << 40) |
						((uint64_t)(ubuf[6]) << 32) |
						((uint64_t)(ubuf[5]) << 24) |
						((uint64_t)(ubuf[4]) << 16) |
						((uint64_t)(ubuf[3]) << 8) | ubuf[2];
		}
		else
		{
			buf += BYTES_PER_BLOCK;
			continue;
		}
		if (frequency == (uint64_t)(FREQ_ONE_MHZ * frequencies[0]))
		{
			if (sweep_started)
			{				
				sweep_count++;

				if (timestamp_normalized == true)
				{
					// set the timestamp of the next sweep
					gettimeofday(&usb_transfer_time, NULL);
				}

				if (one_shot)
				{
					do_exit = true;
				}
				else if (finite_mode && sweep_count == num_sweeps)
				{
					do_exit = true;
				}
			}
			sweep_started = true;
		}
		if (do_exit)
		{
			return 0;
		}
		if (!sweep_started)
		{
			buf += BYTES_PER_BLOCK;
			continue;
		}
		if ((FREQ_MAX_MHZ * FREQ_ONE_MHZ) < frequency)
		{
			buf += BYTES_PER_BLOCK;
			continue;
		}
		/* copy to fftwIn as floats */
		buf += BYTES_PER_BLOCK - (fftSize * 2);
		
		for (i = 0; i < fftSize; i++)
		{
			fftwIn[i][0] = buf[i * 2] * window[i] * 1.0f / 128.0f;
			fftwIn[i][1] = buf[i * 2 + 1] * window[i] * 1.0f / 128.0f;
		}

		buf += fftSize * 2;
		fftwf_execute(fftwPlan);

		for (i = 0; i < fftSize; i++)
		{
			pwr[i] = logPower(fftwOut[i], 1.0f / fftSize);
		}
		
		time_t time_stamp_seconds = usb_transfer_time.tv_sec;
		fft_time = localtime(&time_stamp_seconds);
		strftime(time_str, 50, "%Y-%m-%d, %H:%M:%S", fft_time);
		
		pthread_mutex_lock(&mutex);
		
		float alpha = 2.0/(average+1.0);

		for (i = 0; i < fftSize; i++)
		{
			// Disabled outfile
			//fprintf(outfile, ", %.2f", pwr[i]);			
			int frequency_array_bin = (frequency / 6000) + i;

			if (first_frequency_array_bin == 0)
				first_frequency_array_bin = frequency_array_bin;

			if (average < 1)
			{
				baseline[frequency_array_bin] = pwr[i];
			}
			else
			{				
				baseline[frequency_array_bin] = alpha * pwr[i] + (1 - alpha) * baseline[frequency_array_bin];
				
				if (i == 0 && first_frequency_array_bin == frequency_array_bin)
					average_count++;
				
				if (i == 0 && average_count > average)
				{
					fprintf(stderr, "Averaged %d values\n", average);
					average_count = 0;
				}
			}

			if (baseline_saved)
			{
				float diff = fabs(baseline[frequency_array_bin] - saved_baseline[frequency_array_bin]);
				
				if (diff > baseline_alert_offset)
				{
					float freq_hz = (float)(frequency-(fftSize/2+i));
					float freq_mhz = (float)(frequency-(fftSize/2+i))/1000000;
					int freq_index = -1;

					if (p_records != 0)
						freq_index = lookup_record(p_records, freq_hz);

					fprintf(stderr, "Alert [%f] MHz, power [%f], FCC entry:[%s - %s], baseline differnce %f, sweep count: %u\n",
					freq_mhz, 
					baseline[frequency_array_bin], 
					freq_index == -1 ? "not found" : p_records[freq_index].service,
					freq_index == -1 ? "not found" : p_records[freq_index].notes,
					 diff,
					 sweep_count);
				}
			}
			else if (baseline[frequency_array_bin] > threshold)
			{
				float freq_mhz = (float)(frequency-(fftSize/2+i))/1000000;

				fprintf(stderr, "Threshold alert at freq %f MHz, power %f, threshold %d, sweep count: %u\n",
				 freq_mhz, baseline[frequency_array_bin], threshold, sweep_count);
			}
		}

		pthread_mutex_unlock(&mutex);	
	}
	return 0;
}

static void usage()
{
	fprintf(stderr,
			"Usage:\n"
			"\t[-h] # this help\n"
			"\t[-d serial_number] # Serial number of desired HackRF\n"
			"\t[-t threshold] # Threshold in dB to trigger and alert. Default is -55 dB\n"
			"\t[-a amp_enable] # RX RF amplifier 1=Enable, 0=Disable\n"
			"\t[-f freq_min:freq_max] # minimum and maximum frequencies in MHz\n"
			"\t[-p antenna_enable] # Antenna port power, 1=Enable, 0=Disable\n"
			"\t[-l gain_db] # RX LNA (IF) gain, 0-40dB, 8dB steps\n"
			"\t[-g gain_db] # RX VGA (baseband) gain, 0-62dB, 2dB steps\n"
			"\t[-w bin_width] # FFT bin width (frequency resolution) in Hz, 2445-5000000\n"
			"\t[-W wisdom_file] # Use FFTW wisdom file (will be created if necessary)\n"
			"\t[-P estimate|measure|patient|exhaustive] # FFTW plan type, default is 'measure'\n"
			"\t[-1] # one shot mode\n"
			"\t[-N num_sweeps] # Number of sweeps to perform\n"
			"\t[-B] # binary output\n"
			"\t[-n] # keep the same timestamp within a sweep\n"
			"\t-r filename # output file\n"
			"\n"
			"Output fields:\n"
			"\tdate, time, hz_low, hz_high, hz_bin_width, num_samples, dB, dB, . . .\n");
}

#ifdef _MSC_VER
BOOL WINAPI sighandler(int signum)
{
	if (CTRL_C_EVENT == signum)
	{
		fprintf(stderr, "Caught signal %d\n", signum);
		do_exit = true;
		return TRUE;
	}
	return FALSE;
}
#else
void sigint_callback_handler(int signum)
{
	fprintf(stderr, "Caught signal %d\n", signum);
	do_exit = true;
}
#endif

int import_wisdom(const char *path)
{
	// Returns nonzero
	if (!fftwf_import_wisdom_from_filename(path))
	{
		fprintf(stderr,
				"Wisdom file %s not found; will attempt to create it\n",
				path);
		return 0;
	}

	return 1;
}

int import_default_wisdom()
{
	return fftwf_import_system_wisdom();
}

int export_wisdom(const char *path)
{
	if (path != NULL)
	{
		if (!fftwf_export_wisdom_to_filename(path))
		{
			fprintf(stderr, "Could not write FFTW wisdom file to %s", path);
			return 0;
		}
	}

	return 1;
}



int main(int argc, char **argv)
{
	int opt, i, result = 0;
	const char *path = NULL;
	const char *serial_number = NULL;
	int exit_code = EXIT_SUCCESS;
	struct timeval time_now;
	struct timeval time_prev;
	float time_diff;
	float sweep_rate = 0;
	unsigned int lna_gain = 16, vga_gain = 20;
	uint32_t freq_min = 0;
	uint32_t freq_max = 6000;
	uint32_t requested_fft_bin_width;
	const char *fftwWisdomPath = NULL;
	int fftw_plan_type = FFTW_MEASURE;	
	p_records = read_table();
	
	int index = 0;
	
	for (index = 0; index < BASELINE_SIZE; index++)
	{
		baseline[index] = -80.0;
	}

	while ((opt = getopt(argc, argv, "a:f:p:l:g:d:N:w:W:P:n1BIr:t:h?")) != EOF)
	{
		result = HACKRF_SUCCESS;
		switch (opt)
		{
		case 'd':
			serial_number = optarg;
			break;

		case 'a':
			amp = true;
			result = parse_u32(optarg, &amp_enable);
			break;

		case 'f':
			result = parse_u32_range(optarg, &freq_min, &freq_max);

			if (freq_min >= freq_max)
			{
				fprintf(stderr,
						"argument error: freq_max must be greater than freq_min.\n");
				usage();
				return EXIT_FAILURE;
			}
			if (FREQ_MAX_MHZ < freq_max)
			{
				fprintf(stderr,
						"argument error: freq_max may not be higher than %u.\n",
						FREQ_MAX_MHZ);
				usage();
				return EXIT_FAILURE;
			}
			if (MAX_SWEEP_RANGES <= num_ranges)
			{
				fprintf(stderr,
						"argument error: specify a maximum of %u frequency ranges.\n",
						MAX_SWEEP_RANGES);
				usage();
				return EXIT_FAILURE;
			}

			frequencies[2 * num_ranges] = (uint16_t)freq_min;
			frequencies[2 * num_ranges + 1] = (uint16_t)freq_max;
			num_ranges++;
			break;

		case 'p':
			antenna = true;
			result = parse_u32(optarg, &antenna_enable);
			break;

		case 'l':
			result = parse_u32(optarg, &lna_gain);
			break;

		case 'g':
			result = parse_u32(optarg, &vga_gain);
			break;

		case 'N':
			finite_mode = true;
			result = parse_u32(optarg, &num_sweeps);
			break;

		case 'w':
			result = parse_u32(optarg, &requested_fft_bin_width);
			fftSize = DEFAULT_SAMPLE_RATE_HZ / requested_fft_bin_width;
			break;

		case 'W':
			fftwWisdomPath = optarg;
			break;

		case 'P':
			if (strcmp("estimate", optarg) == 0)
			{
				fftw_plan_type = FFTW_ESTIMATE;
			}
			else if (strcmp("measure", optarg) == 0)
			{
				fftw_plan_type = FFTW_MEASURE;
			}
			else if (strcmp("patient", optarg) == 0)
			{
				fftw_plan_type = FFTW_PATIENT;
			}
			else if (strcmp("exhaustive", optarg) == 0)
			{
				fftw_plan_type = FFTW_EXHAUSTIVE;
			}
			else
			{
				fprintf(stderr, "Unknown FFTW plan type '%s'\n", optarg);
				return EXIT_FAILURE;
			}
			break;

		case 'n':
			timestamp_normalized = true;
			break;

		case '1':
			one_shot = true;
			break;

		case 'B':
			binary_output = true;
			break;

		case 'r':
			path = optarg;
			break;

		case 't':
			result = parse_int(optarg, &threshold);			
			break;

		case 'h':
		case '?':
			usage();
			return EXIT_SUCCESS;

		default:
			fprintf(stderr, "unknown argument '-%c %s'\n", opt, optarg);
			usage();
			fprintf(stderr,
					"argument error: '-%c %s' %s (%d)\n",
					opt,
					optarg,
					hackrf_error_name(result),
					result);
			usage();
			return EXIT_FAILURE;
		}
	}

	// Try to load a wisdom file if specified, otherwise
	// try to load the system-wide wisdom file
	if (fftwWisdomPath)
	{
		import_wisdom(fftwWisdomPath);
	}
	else
	{
		import_default_wisdom();
	}

	if (lna_gain % 8)
	{
		fprintf(stderr, "warning: lna_gain (-l) must be a multiple of 8\n");
	}

	if (vga_gain % 2)
	{
		fprintf(stderr, "warning: vga_gain (-g) must be a multiple of 2\n");
	}

	if (amp)
	{
		if (amp_enable > 1)
		{
			fprintf(stderr, "argument error: amp_enable shall be 0 or 1.\n");
			usage();
			return EXIT_FAILURE;
		}
	}

	if (antenna)
	{
		if (antenna_enable > 1)
		{
			fprintf(stderr,
					"argument error: antenna_enable shall be 0 or 1.\n");
			usage();
			return EXIT_FAILURE;
		}
	}

	if (0 == num_ranges)
	{
		frequencies[0] = (uint16_t)freq_min;
		frequencies[1] = (uint16_t)freq_max;
		num_ranges++;
	}

	/*
	 * The FFT bin width must be no more than a quarter of the sample rate
	 * for interleaved mode. With our fixed sample rate of 20 Msps, that
	 * results in a maximum bin width of 5000000 Hz.
	 */
	if (4 > fftSize)
	{
		fprintf(stderr,
				"argument error: FFT bin width (-w) must be no more than 5000000\n");
		return EXIT_FAILURE;
	}

	/*
	 * The maximum number of FFT bins we support is equal to the number of
	 * samples in a block. Each block consists of 16384 bytes minus 10
	 * bytes for the frequency header, leaving room for 8187 two-byte
	 * samples. As we pad fftSize up to the next odd multiple of four, this
	 * makes our maximum supported fftSize 8180.  With our fixed sample
	 * rate of 20 Msps, that results in a minimum bin width of 2445 Hz.
	 */
	if (8180 < fftSize)
	{
		fprintf(stderr,
				"argument error: FFT bin width (-w) must be no less than 2445\n");
		return EXIT_FAILURE;
	}

	/* In interleaved mode, the FFT bin selection works best if the total
	 * number of FFT bins is equal to an odd multiple of four.
	 * (e.g. 4, 12, 20, 28, 36, . . .)
	 */
	while ((fftSize + 4) % 8)
	{
		fftSize++;
	}

	fft_bin_width = (double)DEFAULT_SAMPLE_RATE_HZ / fftSize;
	fftwIn = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * fftSize);
	fftwOut = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * fftSize);
	fftwPlan = fftwf_plan_dft_1d(fftSize, fftwIn, fftwOut, FFTW_FORWARD, fftw_plan_type);
	pwr = (float *)fftwf_malloc(sizeof(float) * fftSize);
	window = (float *)fftwf_malloc(sizeof(float) * fftSize);

	for (i = 0; i < fftSize; i++)
	{
		window[i] = (float)(0.5f * (1.0f - cos(2 * M_PI * i / (fftSize - 1))));
	}

	/* Execute the plan once to make sure it's ready to go when real
	 * data starts to flow.  See issue #1366
	 */
	fftwf_execute(fftwPlan);

	// reset the timestamp
	memset(&usb_transfer_time, 0, sizeof(usb_transfer_time));

#ifdef _MSC_VER
	if (binary_output)
	{
		_setmode(_fileno(stdout), _O_BINARY);
	}
#endif

	result = hackrf_init();

	if (result != HACKRF_SUCCESS)
	{
		fprintf(stderr,
				"hackrf_init() failed: %s (%d)\n",
				hackrf_error_name(result),
				result);
		usage();
		return EXIT_FAILURE;
	}

	result = hackrf_open_by_serial(serial_number, &device);

	if (result != HACKRF_SUCCESS)
	{
		fprintf(stderr,
				"hackrf_open() failed: %s (%d)\n",
				hackrf_error_name(result),
				result);
		usage();
		return EXIT_FAILURE;
	}

	if ((NULL == path) || (strcmp(path, "-") == 0))
	{
		outfile = stdout;
	}
	else
	{
		outfile = fopen(path, "wb");
	}

	if (NULL == outfile)
	{
		fprintf(stderr, "Failed to open file: %s\n", path);
		return EXIT_FAILURE;
	}
	/* Change outfile buffer to have bigger one to store or read data on/to HDD */
	result = setvbuf(outfile, NULL, _IOFBF, FD_BUFFER_SIZE);
	if (result != 0)
	{
		fprintf(stderr, "setvbuf() failed: %d\n", result);
		usage();
		return EXIT_FAILURE;
	}

#ifdef _MSC_VER
	SetConsoleCtrlHandler((PHANDLER_ROUTINE)sighandler, TRUE);
#else
	signal(SIGINT, &sigint_callback_handler);
	signal(SIGILL, &sigint_callback_handler);
	signal(SIGFPE, &sigint_callback_handler);
	signal(SIGSEGV, &sigint_callback_handler);
	signal(SIGTERM, &sigint_callback_handler);
	signal(SIGABRT, &sigint_callback_handler);
#endif
	fprintf(stderr,
			"call hackrf_sample_rate_set(%.03f MHz)\n",
			((float)DEFAULT_SAMPLE_RATE_HZ / (float)FREQ_ONE_MHZ));
	result = hackrf_set_sample_rate_manual(device, DEFAULT_SAMPLE_RATE_HZ, 1);

	if (result != HACKRF_SUCCESS)
	{
		fprintf(stderr,
				"hackrf_sample_rate_set() failed: %s (%d)\n",
				hackrf_error_name(result),
				result);
		usage();
		return EXIT_FAILURE;
	}

	fprintf(stderr,
			"call hackrf_baseband_filter_bandwidth_set(%.03f MHz)\n",
			((float)DEFAULT_BASEBAND_FILTER_BANDWIDTH / (float)FREQ_ONE_MHZ));
	result = hackrf_set_baseband_filter_bandwidth(
		device,
		DEFAULT_BASEBAND_FILTER_BANDWIDTH);

	if (result != HACKRF_SUCCESS)
	{
		fprintf(stderr,
				"hackrf_baseband_filter_bandwidth_set() failed: %s (%d)\n",
				hackrf_error_name(result),
				result);
		usage();
		return EXIT_FAILURE;
	}

	result = hackrf_set_vga_gain(device, vga_gain);
	result |= hackrf_set_lna_gain(device, lna_gain);

	/*
	 * For each range, plan a whole number of tuning steps of a certain
	 * bandwidth. Increase high end of range if necessary to accommodate a
	 * whole number of steps, minimum 1.
	 */
	for (i = 0; i < num_ranges; i++)
	{
		step_count =
			1 + (frequencies[2 * i + 1] - frequencies[2 * i] - 1) / TUNE_STEP;
		frequencies[2 * i + 1] =
			(uint16_t)(frequencies[2 * i] + step_count * TUNE_STEP);
		fprintf(stderr,
				"Sweeping from %u MHz to %u MHz\n",
				frequencies[2 * i],
				frequencies[2 * i + 1]);
	}

	result = hackrf_init_sweep(
		device,
		frequencies,
		num_ranges,
		BYTES_PER_BLOCK,
		TUNE_STEP * FREQ_ONE_MHZ,
		OFFSET,
		INTERLEAVED);

	if (result != HACKRF_SUCCESS)
	{
		fprintf(stderr,
				"hackrf_init_sweep() failed: %s (%d)\n",
				hackrf_error_name(result),
				result);
		return EXIT_FAILURE;
	}

	result |= hackrf_start_rx_sweep(device, rx_callback, NULL);

	if (result != HACKRF_SUCCESS)
	{
		fprintf(stderr,
				"hackrf_start_rx_sweep() failed: %s (%d)\n",
				hackrf_error_name(result),
				result);
		usage();
		return EXIT_FAILURE;
	}

	if (amp)
	{
		fprintf(stderr, "call hackrf_set_amp_enable(%u)\n", amp_enable);
		result = hackrf_set_amp_enable(device, (uint8_t)amp_enable);
		if (result != HACKRF_SUCCESS)
		{
			fprintf(stderr,
					"hackrf_set_amp_enable() failed: %s (%d)\n",
					hackrf_error_name(result),
					result);
			usage();
			return EXIT_FAILURE;
		}
	}

	if (antenna)
	{
		fprintf(stderr, "call hackrf_set_antenna_enable(%u)\n", antenna_enable);
		result = hackrf_set_antenna_enable(device, (uint8_t)antenna_enable);
		if (result != HACKRF_SUCCESS)
		{
			fprintf(stderr,
					"hackrf_set_antenna_enable() failed: %s (%d)\n",
					hackrf_error_name(result),
					result);
			usage();
			return EXIT_FAILURE;
		}
	}

	gettimeofday(&t_start, NULL);
	time_prev = t_start;

	fprintf(stderr, "Stop with Ctrl-C\n");

	while ((hackrf_is_streaming(device) == HACKRF_TRUE) && (do_exit == false))
	{
		float time_difference;
		m_sleep(50);

		if (kbhit())
		{
			int key_pressed = getchar();
			fprintf(stderr, "\nKey pressed: %c.\n", key_pressed);

			process_command(key_pressed);
		}

		gettimeofday(&time_now, NULL);

		if (TimevalDiff(&time_now, &time_prev) >= 1.0f)
		{
			time_difference = TimevalDiff(&time_now, &t_start);
			sweep_rate = (float)sweep_count / time_difference;
			
			// fprintf(stderr,					
			// 		"%u total sweeps completed, %.2f sweeps/second\n",
			// 		sweep_count,
			// 		sweep_rate);

			if (byte_count == 0)
			{
				exit_code = EXIT_FAILURE;
				fprintf(stderr,
						"\nCouldn't transfer any data for one second.\n");
				break;
			}
			byte_count = 0;
			time_prev = time_now;
		}
	}

	fflush(outfile);
	result = hackrf_is_streaming(device);

	if (do_exit)
	{
		fprintf(stderr, "\nExiting...\n");

		if (p_records != 0)
			free (p_records);
	}
	else
	{
		fprintf(stderr,
				"\nExiting... hackrf_is_streaming() result: %s (%d)\n",
				hackrf_error_name(result),
				result);
	}

	gettimeofday(&time_now, NULL);
	time_diff = TimevalDiff(&time_now, &t_start);
	
	if ((sweep_rate == 0) && (time_diff > 0))
	{
		sweep_rate = sweep_count / time_diff;
	}
	
	fprintf(stderr,
			"Total sweeps: %u in %.5f seconds (%.2f sweeps/second)\n",
			sweep_count,
			time_diff,
			sweep_rate);

	if (device != NULL)
	{
		result = hackrf_close(device);
		if (result != HACKRF_SUCCESS)
		{
			fprintf(stderr,
					"hackrf_close() failed: %s (%d)\n",
					hackrf_error_name(result),
					result);
		}
		else
		{
			fprintf(stderr, "hackrf_close() done\n");
		}

		hackrf_exit();
		fprintf(stderr, "hackrf_exit() done\n");
	}

	fflush(outfile);
	if ((outfile != NULL) && (outfile != stdout))
	{
		fclose(outfile);
		outfile = NULL;
		fprintf(stderr, "fclose() done\n");
	}
	fftwf_free(fftwIn);
	fftwf_free(fftwOut);
	fftwf_free(pwr);
	fftwf_free(window);
	fftwf_free(ifftwIn);
	fftwf_free(ifftwOut);
	export_wisdom(fftwWisdomPath);
	fprintf(stderr, "exit\n");
	return exit_code;
}
