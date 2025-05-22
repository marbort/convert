import requests

def search_stop_places():
    url = "https://api.entur.org/stop_places/1.0/graphql/"
    query =""" 
    {    
         topographicPlace(query: "Nesbru") {
        id
        name {
            value
        }
    }
}
"""
    headers = {
        "Content-Type": "application/json",
        "ET-Client-Name": "your-client-id"  # Use your own identifier
    }
    response = requests.get(url,json={"query": query}, headers=headers)
    data=response.json()
    print(data)

search_stop_places()