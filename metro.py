import requests
from datetime import date, datetime, timedelta
from colorama import Fore, Style, Back



# Define the GraphQL query
stops={"Blindern":6332,"Kalbakken":5810}
line_colors={"1":Back.CYAN,"2":Back.RED,"3":Back.MAGENTA,"4":Back.BLUE,"5":Back.GREEN}
query = """
{
  stopPlace(id: "NSR:StopPlace:##STOP_ID##") {
    name
    estimatedCalls(timeRange: 7200, numberOfDepartures: 10) {
      realtime
      aimedDepartureTime
      expectedDepartureTime
      destinationDisplay {
        frontText
      }
      serviceJourney {
        line {
          publicCode
          transportMode
        }
      }
    }
  }
}
"""

# Set the endpoint and headers
url = "https://api.entur.io/journey-planner/v3/graphql"
headers = {
    "Content-Type": "application/json",
    "ET-Client-Name": "your-client-id"  # Use your own identifier
}

# Send the request
responses = {}
for i in stops:
    responses[i] = requests.post(url, json={"query": query.replace('##STOP_ID##',str(stops[i]))}, headers=headers)
#print(responses)
# Print result
for i in responses:
    if responses[i].status_code == 200:
        data = responses[i].json()
        #print(data)
        print(f"DEPARTURES FROM {i}")
        for call in data["data"]["stopPlace"]["estimatedCalls"]:
            date_format = '%Y-%m-%dT%H:%M:%S%z'
            time=datetime.strptime(call["aimedDepartureTime"],date_format)
            time_expected=datetime.strptime(call["expectedDepartureTime"],date_format)
            delay=time_expected-time
            
            print(f"{line_colors[call['serviceJourney']['line']['publicCode']]}Line {call['serviceJourney']['line']['publicCode']:<4} to {call['destinationDisplay']['frontText']:<30} at {time.strftime('%H:%M')} expected at {time_expected.strftime('%H:%M')} delay {int(delay.total_seconds()//60):>3} minutes{Style.RESET_ALL}")
            
        
        
        #print(f"Line {call['serviceJourney']['line']['publicCode']} to {call['destinationDisplay']['frontText']} at {call['expectedDepartureTime']}")
        
        #print(f"Line {call['serviceJourney']['line']['publicCode']} to {call['destinationDisplay']['frontText']} at {call['aimedDepartureTime']} expected at  {call['expectedDepartureTime']}")
    else:
        print("FError:", responses[i].status_code)

