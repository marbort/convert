# Import Module
from tkinter import *
import requests
from datetime import date, datetime, timedelta
from colorama import Fore, Style, Back
from tkinter import Tk, font


def clicked():
    responses,delay,departures,stops=metro_times()
    update_times(responses,delay,departures,stops)

def update_times(responses,delay,departures,stops):
    for j in range(len(stops)):
        for i in range(len(label_sizes)):
            for departure in range(5):
                lbl = Label(root, text = departures[list(stops.keys())[j]][departure][slices1[i]:slices2[i]], 
                            font=("Helvetica", 22, 'bold'),justify='left',width=label_sizes[i],anchor='w',
                            fg=colors[departures[list(stops.keys())[j]][departure][:6]])
                lbl.grid(column=i, row=9*j+3+departure,sticky='ew')
                root.grid_rowconfigure(9*j+1, weight=1)
        root.grid_columnconfigure(i, weight=1)

def metro_times():
    # Define the GraphQL query
    stops={"Blindern":6332,"Kalbakken":5810,"Oslo S":3990}
    line_colors={"1":Fore.CYAN,"2":Fore.RED,"3":Fore.MAGENTA,"4":Fore.BLUE,"5":Fore.GREEN}
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
    departures = {}
    for i in stops:
        responses[i] = requests.post(url, json={"query": query.replace('##STOP_ID##',str(stops[i]))}, headers=headers)
        departures[i]=[]
    # Print result
    for i,name in enumerate(responses):
        if responses[name].status_code == 200:
            data = responses[name].json()
            for call in data["data"]["stopPlace"]["estimatedCalls"]:
                date_format = '%Y-%m-%dT%H:%M:%S%z'
                time=datetime.strptime(call["aimedDepartureTime"],date_format)
                time_expected=datetime.strptime(call["expectedDepartureTime"],date_format)
                delay=time_expected-time
                departures[list(stops.keys())[i]].append(f"Line {call['serviceJourney']['line']['publicCode']:<4}  {call['destinationDisplay']['frontText']:<30} at {time.strftime('%H:%M')} expected at {time_expected.strftime('%H:%M')} delay {int(delay.total_seconds()//60):>3} minutes")
                
            
            
            #print(f"Line {call['serviceJourney']['line']['publicCode']} to {call['destinationDisplay']['frontText']} at {call['expectedDepartureTime']}")
            
            #print(f"Line {call['serviceJourney']['line']['publicCode']} to {call['destinationDisplay']['frontText']} at {call['aimedDepartureTime']} expected at  {call['expectedDepartureTime']}")
        else:
            print("FError:", responses[i].status_code)
    return(responses,str(delay),departures,stops)


label_sizes=[35,15,15,15]
slices1=[0,45,63,77]
slices2=[42,50,68,87]
label_headers=["Line","Time","Expected","Delay"]
colors={"Line 1":"cyan","Line 2":"orange","Line 3":"purple","Line 4":"blue","Line 5":"green"}
#f = tkFont.Font(family='helvetica', size=20)

# create root window
root = Tk()

# root window title and dimension
root.title("Metro Live Times")
# Set geometry(widthxheight)
root.geometry('1920x1080')
responses,delay,departures,stops=metro_times()
print(departures)

for j in range(len(stops)):
    lbl = Label(root, text =f"Departures from {list(stops.keys())[j]}", justify='left', font=("Helvetica", 24))
    lbl.grid(column=0, row=9*j,columnspan=len(label_sizes),sticky='W')
    for i in range(len(label_sizes)):
        lbl = Label(root, text = label_headers[i], font=("Helvetica", 22),justify='left',width=label_sizes[i],anchor='w')
        lbl.grid(column=i, row=9*j+1,sticky='ew')
        for departure in range(5):
                lbl = Label(root, text = departures[list(stops.keys())[j]][departure][slices1[i]:slices2[i]], 
                            font=("Helvetica", 22, 'bold'),justify='left',width=label_sizes[i],anchor='w',
                            fg=colors[departures[list(stops.keys())[j]][departure][:6]])
                lbl.grid(column=i, row=9*j+3+departure,sticky='ew')
                root.grid_rowconfigure(9*j+3+departure, weight=1)
        root.grid_columnconfigure(i, weight=1)
        
    #CODE FOR TIMES
    lbl = Label(root, text ="----------------", justify='left', font=("Helvetica", 20))
    lbl.grid(column=0, row=9*j+8,columnspan=len(label_sizes),sticky='W')

print(root.grid_size())
# button widget with red color text
# inside
btn = Button(root, text = "Update Times" ,
             fg = "black", command=clicked,  font=("Helvetica", 20))
# set Button grid
btn.grid(column=0, row=root.grid_size()[-1]+1,columnspan=len(label_sizes),sticky='ew')


#for i in range(len(responses)):
    




# Execute Tkinter
root.mainloop()