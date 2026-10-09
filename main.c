// ============================================================================
// main.c -- Interactive DNA Sequence Analyzer for PicoRV32 SoC
//
// Operations Supported:
//   Mode 1: Sequence Matching (Whole-sequence exact equality comparison)
//   Mode 2: Motif Detection   (Exact substring search across candidate positions)
//
// Memory Map:
//   0x10000000 : Application Base (Volatile Block RAM, 8 KB)
//   0x20000000 : UART DATA (read: RX byte, write: TX byte)
//   0x20000004 : UART STATUS (bit 0 = tx_busy, bit 1 = rx_valid)
//   0x30000000 : LED Register (16-bit)
//   0x40000000 : DNA MMIO Base
//     +0x00 : REF_DATA       (W)   - Write 16 packed bases into ref_mem
//     +0x04 : MOTIF_DATA     (W)   - Write 16 packed bases into motif_mem
//     +0x08 : CONFIG         (W/R) - [7:0]=ref_len, [15:8]=motif_len
//     +0x0C : CONTROL        (W/R) - bit 0=START, bit 1=RESET_PTRS, bit 2=MODE
//     +0x10 : STATUS         (R)   - bit 0=BUSY, bit 1=DONE
//     +0x14 : MATCH_COUNT    (R)   - [7:0]=Matches found (or 1=FOUND/0=NOT FOUND)
//     +0x18 : POSITION_INDEX (W/R) - [6:0]=Desired match index to read
//     +0x1C : POSITION_DATA  (R)   - [7:0]=Match starting position
//
// DNA Encoding:
//   A = 2'b00, C = 2'b01, G = 2'b10, T = 2'b11
// ============================================================================

#define UART_DATA       (*(volatile unsigned int *)0x20000000)
#define UART_STATUS     (*(volatile unsigned int *)0x20000004)
#define LED_REG         (*(volatile unsigned int *)0x30000000)

#define DNA_REF_DATA    (*(volatile unsigned int *)0x40000000)
#define DNA_MOTIF_DATA  (*(volatile unsigned int *)0x40000004)
#define DNA_CONFIG      (*(volatile unsigned int *)0x40000008)
#define DNA_CONTROL     (*(volatile unsigned int *)0x4000000C)
#define DNA_STATUS      (*(volatile unsigned int *)0x40000010)
#define DNA_MATCH_COUNT (*(volatile unsigned int *)0x40000014)
#define DNA_POS_INDEX   (*(volatile unsigned int *)0x40000018)
#define DNA_POS_DATA    (*(volatile unsigned int *)0x4000001C)

#define UART_TX_BUSY    0x1
#define UART_RX_VALID   0x2

#define DNA_STATUS_BUSY 0x1
#define DNA_STATUS_DONE 0x2

#define DNA_START       0x1
#define DNA_RESET_PTRS  0x2
#define DNA_MODE_SEQ    (0 << 2)  // Mode 0: Sequence Matching
#define DNA_MODE_MOTIF  (1 << 2)  // Mode 1: Motif Detection

// LED Status Patterns (preserving bootloader meanings 0x0001, 0x0003, 0xFFFF)
#define LED_APP_MENU       0x0010  // Main menu / idle
#define LED_WAIT_TARGET    0x0020  // Waiting for target sequence
#define LED_WAIT_QUERY     0x0040  // Waiting for query/motif sequence
#define LED_PROCESSING     0x0080  // Hardware accelerator running
#define LED_SEQ_FOUND      0x0101  // Mode 1: Sequence Match FOUND
#define LED_SEQ_NOT_FOUND  0x0100  // Mode 1: Sequence Match NOT FOUND
#define LED_MOTIF_COMPLETE 0x0200  // Mode 2: Motif complete (+ count in lower 8 bits)



static char uart_getc(void)
{
    while (!(UART_STATUS & UART_RX_VALID))
        ;
    return (char)(UART_DATA & 0xFF);
}

static void uart_putc(char c)
{
    while (UART_STATUS & UART_TX_BUSY)
        ;
    UART_DATA = (unsigned char)c;
}

static void uart_puts(const char *s)
{
    while (*s) {
        if (*s == '\n') {
            while (UART_STATUS & UART_TX_BUSY)
                ;
            UART_DATA = '\r';
        }
        while (UART_STATUS & UART_TX_BUSY)
            ;
        UART_DATA = (unsigned char)*s++;
    }
}

// Decimal printing without software division library dependency
static void uart_put_uint(unsigned int n)
{
    if (n == 0) {
        uart_putc('0');
        return;
    }

    char buf[12];
    int idx = 0;

    while (n > 0) {
        unsigned int q = 0;
        unsigned int r = n;
        while (r >= 10) {
            r -= 10;
            q++;
        }
        buf[idx++] = (char)('0' + r);
        n = q;
    }

    for (int i = idx - 1; i >= 0; i--) {
        uart_putc(buf[i]);
    }
}


static int read_dna_sequence(char *buf, int max_len)
{
    int len = 0;

    while (1) {
        char c = uart_getc();

        // Handle Carriage Return / Line Feed (Enter key)
        if (c == '\r' || c == '\n') {
            if (len > 0) {
                uart_puts("\r\n");
                buf[len] = '\0';
                return len;
            }
            continue; // Ignore empty lines
        }

        // Handle Backspace (ASCII 8 or DEL 127)
        if (c == '\b' || c == 127) {
            if (len > 0) {
                len--;
                uart_puts("\b \b"); // Erase character from terminal display
            }
            continue;
        }

        // Convert lowercase to uppercase
        if (c >= 'a' && c <= 'z') {
            c -= ('a' - 'A');
        }

        // Validate DNA character: A, C, G, T
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            if (len < max_len) {
                buf[len++] = c;
                uart_putc(c); // Echo character to terminal
            }
        }
        // Discard invalid characters silently
    }
}


static char read_choice(void)
{
    while (1) {
        char c = uart_getc();
        if (c == '1' || c == '2') {
            uart_putc(c);
            uart_puts("\r\n\n");
            return c;
        }
    }
}


static inline unsigned int encode_base(char c)
{
    switch (c) {
        case 'A': return 0x0;
        case 'C': return 0x1;
        case 'G': return 0x2;
        case 'T': return 0x3;
        default:  return 0x0;
    }
}


static void run_dna_accelerator(int mode, const char *target, int target_len,
                                const char *query, int query_len,
                                unsigned int *out_match_count)
{
    // 1. Reset memory write pointers in MMIO wrapper & set mode
    unsigned int mode_bit = (mode == 1) ? DNA_MODE_MOTIF : DNA_MODE_SEQ;
    DNA_CONTROL = mode_bit | DNA_RESET_PTRS;

    // 2. Stream packed 16-base words to REF_DATA (target sequence)
    int target_words = (target_len + 15) / 16;
    for (int w = 0; w < target_words; w++) {
        unsigned int word = 0;
        for (int b = 0; b < 16; b++) {
            int base_idx = w * 16 + b;
            unsigned int code = 0;
            if (base_idx < target_len) {
                code = encode_base(target[base_idx]);
            }
            word |= (code << (b * 2));
        }
        DNA_REF_DATA = word;
    }

    // 3. Stream packed 16-base words to MOTIF_DATA (query or motif sequence)
    int query_words = (query_len + 15) / 16;
    for (int w = 0; w < query_words; w++) {
        unsigned int word = 0;
        for (int b = 0; b < 16; b++) {
            int base_idx = w * 16 + b;
            unsigned int code = 0;
            if (base_idx < query_len) {
                code = encode_base(query[base_idx]);
            }
            word |= (code << (b * 2));
        }
        DNA_MOTIF_DATA = word;
    }

    // 4. Configure lengths in CONFIG register:
    //    [7:0] = ref_len, [15:8] = motif_len
    DNA_CONFIG = ((unsigned int)target_len & 0xFF) | (((unsigned int)query_len & 0xFF) << 8);

    // 5. Trigger START pulse with mode selection
    DNA_CONTROL = mode_bit | DNA_START;

    // 6. Wait for hardware accelerator completion (poll STATUS register)
    while (!(DNA_STATUS & DNA_STATUS_DONE))
        ;

    // 7. Read hardware result
    *out_match_count = DNA_MATCH_COUNT & 0xFF;
}


int main(void)
{
    char target_seq[130];
    char query_seq[130];

    LED_REG = LED_APP_MENU;

    while (1) {
        LED_REG = LED_APP_MENU;

        uart_puts("========================================\n");
        uart_puts("       DNA SEQUENCE ANALYZER\n");
        uart_puts("========================================\n\n");
        uart_puts("Select Operation:\n\n");
        uart_puts("1. Sequence Matching\n");
        uart_puts("2. Motif Detection\n\n");
        uart_puts("Enter choice: ");

        char choice = read_choice();

        if (choice == '1') {


            LED_REG = LED_WAIT_TARGET;
            uart_puts("Enter target sequence:\n");
            int target_len = read_dna_sequence(target_seq, 128);

            LED_REG = LED_WAIT_QUERY;
            uart_puts("Enter query sequence:\n");
            int query_len = read_dna_sequence(query_seq, 128);

            LED_REG = LED_PROCESSING;
            uart_puts("Processing...\n\n");

            // Execute hardware comparison on FPGA accelerator
            unsigned int match_result = 0;
            run_dna_accelerator(0, target_seq, target_len, query_seq, query_len, &match_result);

            if (match_result == 1) {
                LED_REG = LED_SEQ_FOUND;
                uart_puts("Sequence Match: FOUND\n\n");
            } else {
                LED_REG = LED_SEQ_NOT_FOUND;
                uart_puts("Sequence Match: NOT FOUND\n\n");
            }

        } else if (choice == '2') {

            LED_REG = LED_WAIT_TARGET;
            uart_puts("Enter target sequence:\n");
            int target_len = read_dna_sequence(target_seq, 128);

            LED_REG = LED_WAIT_QUERY;
            uart_puts("Enter motif sequence:\n");
            int motif_len = read_dna_sequence(query_seq, 16);

            LED_REG = LED_PROCESSING;
            uart_puts("Processing...\n\n");

            // Execute hardware motif search on FPGA accelerator
            unsigned int match_count = 0;
            run_dna_accelerator(1, target_seq, target_len, query_seq, motif_len, &match_count);

            LED_REG = LED_MOTIF_COMPLETE | (match_count & 0xFF);

            uart_puts("Match count: ");
            uart_put_uint(match_count);
            uart_puts("\n");

            uart_puts("Match positions: ");
            if (match_count == 0) {
                uart_puts("None\n\n");
            } else {
                for (unsigned int i = 0; i < match_count; i++) {
                    DNA_POS_INDEX = i;
                    unsigned int pos = DNA_POS_DATA & 0xFF;
                    uart_put_uint(pos);
                    if (i + 1 < match_count) {
                        uart_putc(' ');
                    }
                }
                uart_puts("\n\n");
            }
        }
    }

    return 0;
}
