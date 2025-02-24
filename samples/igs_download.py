import requests
from datetime import datetime, timedelta
from ftplib import FTP, FTP_TLS
import gzip
import os
from urllib.parse import urlparse
import shutil  # Required for file extraction
import subprocess
from tqdm import tqdm  # Progress bar


# Configuration
#
LOCAL_BASE_DIR = "../data"  # Data directory
FILE_LIST = "../data/igs_files.txt"   # File containing destination & FTP URLs


def gps_to_datetime(gps_week, gps_day_of_week):
    # GPS epoch: January 6, 1980
    gps_epoch = datetime(1980, 1, 6)

    # Add weeks and days to the GPS epoch
    date = gps_epoch + timedelta(weeks=gps_week, days=gps_day_of_week)

    return date


def extract_gz(gz_path):
    """Extracts a .gz file to the specified directory."""
    try:
        output_file = os.path.splitext(gz_path)[0]  # Remove .gz extension
        with gzip.open(gz_path, 'rb') as f_in, open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)  # Delete .gz after extraction
        return output_file
    except Exception as e:
        print(f"Warning: Failed to extract {gz_path} - {e}")


def extract_Z(z_path):
    """Extracts a .Z file using the 'gzip' command."""
    try:
        # Rename .Z to .gz
        #
        gz_path = z_path[:-2] + ".gz"
        os.rename(z_path, gz_path)

        # Use the gzip command to extract the file
        #
        subprocess.run(
            ["gzip", "-df", gz_path],  # '-d' to decompress, '-f' to overwrite
            check=True,
            capture_output=True,
            text=True
        )
        return os.path.splitext(gz_path)[0]  # Remove .gz extension
    except FileNotFoundError:
        print("Error: 'gzip' command not found!")
    except Exception as e:
        print(f"Warning: Failed to extract {z_path} - {e}")


def connect_ftp(host):
    """Connect to an FTP server using anonymous login."""
    if 'cddis' in host:
        ftp = FTP_TLS(host)
    else:
        ftp = FTP(host)
    ftp.login()  # Anonymous login
    if 'cddis' in host:
        ftp.prot_p()  # switch to secure data connection
    ftp.voidcmd("TYPE I")  # Switch to binary mode
    print(f"Connected to {host} as anonymous (Binary Mode)")
    return ftp


def read_file_list():
    """Read the list of destination folders and FTP URLs from a file."""
    if not os.path.exists(FILE_LIST):
        print(f"Error: {FILE_LIST} not found!")
        return []

    entries = []
    with open(FILE_LIST, "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if line.startswith('#') or len(parts) == 0:
                continue
            if len(parts) != 2:
                print(f"Skipping invalid line: {line.strip()}")
                continue
            dest_folder, ftp_url = parts
            entries.append((dest_folder, ftp_url))

    print(f"Found {len(entries)} files to download.")
    return entries


def download_ftp(ftp, filename, local_filename):
    """Download file"""
    try:
        total_size = ftp.size(filename)  # Get file size
    except:
        total_size = None  # No progress bar

    if total_size:

        progress_bar = tqdm(total=total_size, unit="B",
                            unit_scale=True, desc=filename)

        def write_chunk(data):
            local_file.write(data)
            progress_bar.update(len(data))

        try:
            with open(local_filename, "wb") as local_file:
                ftp.retrbinary(f"RETR {filename}", write_chunk)
        except Exception as e:
            os.remove(local_filename)
            print(f"Warning: Failed to download {filename} - {e}")

        progress_bar.close()

    else:

        # Download without progress bar if file size is unknown
        try:
            with open(local_filename, "wb") as local_file:
                ftp.retrbinary(f"RETR {filename}", local_file.write)
            print(f"Downloaded: {filename} -> {local_filename}")
        except Exception as e:
            os.remove(local_filename)
            print(f"Warning: Failed to download {filename} - {e}")


def download_http(url, filename, local_filename):
    """Download a file via HTTP/HTTPS."""
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an error for bad status codes

        total_size = int(response.headers.get('content-length', 0))
        progress_bar = tqdm(total=total_size, unit="B",
                            unit_scale=True, desc=filename)

        with open(local_filename, 'wb') as local_file:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    local_file.write(chunk)
                    progress_bar.update(len(chunk))

        progress_bar.close()

    except requests.exceptions.RequestException as e:
        print(f"Warning: Failed to download {filename} from {url} - {e}")


def download_files():
    """Download files from the provided list of FTP URLs."""
    entries = read_file_list()
    if not entries:
        print("No valid entries found. Exiting.")
        return

    current_host = None
    ftp = None

    for dest_folder, url in entries:

        parsed = urlparse(url)
        if parsed.scheme != "ftp" and 'http' not in parsed.scheme:
            print(f"Skipping invalid URL (not FTP/HTTP(S)): {url}")
            continue

        host = parsed.hostname
        remote_path = parsed.path
        filename = os.path.basename(remote_path)
        local_folder = os.path.join(LOCAL_BASE_DIR, dest_folder)
        local_filename = os.path.join(local_folder, filename)

        try:

            # Ensure local folder exists
            #
            os.makedirs(local_folder, exist_ok=True)

            if parsed.scheme == "ftp":

                # Connect to the FTP server if different from the previous one
                #
                if host != current_host:
                    if ftp:
                        ftp.quit()
                    ftp = connect_ftp(host)
                    current_host = host

                # Change to the appropriate directory
                #
                remote_dir = os.path.dirname(remote_path)
                ftp.cwd(remote_dir)

                # Download file
                #
                download_ftp(ftp, filename, local_filename)

            elif "http" in parsed.scheme:
                host = parsed.scheme+"://"+host+remote_path
                download_http(host, filename, local_filename)

            # Check if it is a compressed file and extract it
            #
            if os.path.exists(local_filename):
                if local_filename.lower().endswith(".gz"):
                    local_filename = extract_gz(local_filename)
                elif local_filename.upper().endswith(".Z"):
                    local_filename = extract_Z(local_filename)

            # Short to long RINEX filename
            #
            filename = os.path.basename(local_filename)
            if filename.startswith("COM"):
                epoch = gps_to_datetime(int(filename[3:7]), int(filename[7]))
                if filename.endswith(".BIA"):
                    filename = f"COD0MGXFIN_{epoch:%Y%j0000}_01D_01D_OSB.BIA"
                elif filename.endswith(".CLK"):
                    filename = f"COD0MGXFIN_{epoch:%Y%j0000}_01D_30S_CLK.CLK"
                elif filename.endswith(".EPH"):
                    filename = f"COD0MGXFIN_{epoch:%Y%j0000}_01D_05M_ORB.SP3"
                os.rename(local_filename,
                          os.path.join(local_folder, filename))

        except Exception as e:
            print(f"Warning: Failed to download {filename} from {url} - {e}")

    if ftp:
        ftp.quit()
    print("FTP session closed.")


if __name__ == "__main__":
    download_files()
