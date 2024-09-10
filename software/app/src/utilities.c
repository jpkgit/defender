#include "utilities.h"

// Function to check if a key has been pressed
int kbhit(void)
{
	struct termios oldt, newt;
	int ch;
	int oldf;

	tcgetattr(STDIN_FILENO, &oldt);
	newt = oldt;
	newt.c_lflag &= ~(ICANON | ECHO);
	tcsetattr(STDIN_FILENO, TCSANOW, &newt);
	oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
	fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

	ch = getchar();

	tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
	fcntl(STDIN_FILENO, F_SETFL, oldf);

	if (ch != EOF)
	{
		ungetc(ch, stdin);
		return 1;
	}

	return 0;
}

int parse_int(char *s, int *const value)
{
	int i = atoi(s);
	*value = i;
	return HACKRF_SUCCESS;
}

int parse_u32(char *s, uint32_t *const value)
{
	uint_fast8_t base = 10;
	char *s_end;
	uint64_t ulong_value;

	if (strlen(s) > 2)
	{
		if (s[0] == '0')
		{
			if ((s[1] == 'x') || (s[1] == 'X'))
			{
				base = 16;
				s += 2;
			}
			else if ((s[1] == 'b') || (s[1] == 'B'))
			{
				base = 2;
				s += 2;
			}
		}
	}

	s_end = s;
	ulong_value = strtoul(s, &s_end, base);
	if ((s != s_end) && (*s_end == 0))
	{
		*value = (uint32_t)ulong_value;
		return HACKRF_SUCCESS;
	}
	else
	{
		return HACKRF_ERROR_INVALID_PARAM;
	}
}

int parse_u32_range(char *s, uint32_t *const value_min, uint32_t *const value_max)
{
	int result;

	char *sep = strchr(s, ':');
	if (!sep)
	{
		return HACKRF_ERROR_INVALID_PARAM;
	}

	*sep = 0;

	result = parse_u32(s, value_min);
	if (result != HACKRF_SUCCESS)
	{
		return result;
	}
	result = parse_u32(sep + 1, value_max);
	if (result != HACKRF_SUCCESS)
	{
		return result;
	}

	return HACKRF_SUCCESS;
}