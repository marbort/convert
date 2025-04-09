import os
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import sys

# Function to download all images from a webpage
def download_images_from_webpage(url, save_directory="images"):
    # Create the directory if it doesn't exist
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    # Get the webpage content
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to retrieve the webpage: {response.status_code}")
        return
    
    # Parse the webpage
    soup = BeautifulSoup(response.text, "html.parser")
    
    # Find all image tags
    img_tags = soup.find_all("img")
    
    for img in img_tags:
        img_url = img.get("src")  # Get the source URL of the image
        if not img_url:
            continue
        
        # Handle relative URLs
        img_url = urljoin(url, img_url)
        
        try:
            # Download the image
            img_response = requests.get(img_url, stream=True)
            if img_response.status_code == 200:
                # Extract the image filename
                filename = os.path.join(save_directory, os.path.basename(img_url))
                
                # Save the image to the local directory
                with open(filename, "wb") as file:
                    for chunk in img_response.iter_content(1024):
                        file.write(chunk)
                
                print(f"Downloaded: {filename}")
            else:
                print(f"Failed to download image: {img_url}")
        except Exception as e:
            print(f"Error downloading {img_url}: {e}")

# Example Usage
with open('paths.txt','r') as ifile:
    paths = ifile.readlines()
    print(paths)
    for line in paths:
        download_images_from_webpage(line.strip())
