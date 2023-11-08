
/**** macro variables to be used in new program ****/
/***

datain = original data 
checked_data = input of previous trials that have been checked
eligible = variable for eligible trial 
m = time for trial
n = number of controls for each case 
seed = starting seed for selecting  controls
check_for_cases = indicator if we should stop program at m once the controls have been chosen
recheck_for_cases = indicator if we should check again after adding in new controls. This will allow for an interative process for adding in controls at trial m
    and then checking to see if any were cases.
additional_controls = after controls for trial m have been checked, this is the number of additional controls that need to be selected


data contained in checked_data will be the first m trials that would be used in the final analysis.

When checking the data for new cases, the original data needs to also be fixed to contain any new cases. This will affect trials after the checked trial at time = m. The data listed
in checked data will have the new cases added into trial m and the control removed. This data will be used for adding in the additional trials to preserve those selected in the initial 
selection of controls at trial m. Each additional call to the macro should only add in additional controls on top of those initially selected.
Authors: Roger Logan, Miguel Hernan, Barbra Dickerman, Goodarz Danaei
Written by Roger Logan. rwlogan@hsph.harvard.edu
Version Sept, 2020.

******/
/* for first run we will assume that data has been check for all possible cases */

%macro ccc(datain = , 
           checked_data = , /* collection of trials upto value contained in &m that has been checked for controls that could be cases */
		   id =,
		   time= ,
		   event= ,
		   treatment = ,
           eligible= ,		   
		   eligible_controls = , /* variable or condition for further restricting control selection */
           pp = 1, /* 1 = run pp analysis with weights, 0 = run itt analysis */
		   ccratio= ,

		   check_for_cases = 0, /* only run one iteration of selecting controls for trial contained in m */
		   checked_for_cases = 0, /* data supplied in checked_data was checked and possibly new controls need to be selected */
		   recheck_for_cases = 0, /* used when previously ran with check for cases = 1 to add in new controls and then resubmit the data, needs cheked_for_cases = 1 */
		   additional_controls = ,
		   check_controls = 0, /* when all trials have been created go back and check the covariates for the controls */
		   checked_controls = 0 , /* when the data for the controls have been checked , prior to running final analysis */

           outcomeCov = ,     /* function of baseline covariates used in final model.   */
           outcomeClass = ,   /* Any categorical variables used in the final model */
           outcomeCov_var = , /* list of individual baseline variables used in final model */
           model_var      = , /************************************************************************** 
                                           variables of interest to be used in final model. These variables need to be defined in the sub
                                           macro %create_final_variables which will use the submacro %building_blocks. These two macros are defined
                                           defined for the following two cases when model_var is left to be missing (by defualt)
                                          
                                         
                           
                                         ***************************************************************************/
		   include_expansion_time = 1, /* include for_period and for_period2 in switch models */
           include_expansion_time_outcome = 1, /* include for period in model for outcome */
           include_followup_time = 1 , /* include follow up time in switch model (based on for_period ) */
           include_followup_time_outcome = 1 , /* include follow up time in model for outcome */

           pool_models = 0 , /* pool treatment models over previous value of treatment  */
           condition0= ,
           condition1= ,
           where_var =  ,
		   calculate_var = 0,
		   use_stabilized_weights = 0 ,
		   eligible_wts_0 = ,
		   eligible_wts_1= ,
		   cov_switchd = ,
		   class_switchd = ,
		   cov_switchn = ,
		   class_switchn = ,
		   cense = ,
		   cov_censed = ,
		   class_censed= ,
		   cov_censen = ,
		   class_censen = ,
		   final_analysis = 1,
           run_p99_analysis = 0 ,     /* run the final model with truncating the weights at the 1st and 99th percentile */
                     run_user_limits_analysis = 0, /* run the final model with truncating the weights using user defined limits */                            
                     lower_weight = ,     /* use lower weight instead of the p01 given by the proc univariate, used when p99 = 1 */
                     upper_weight = ,     /* use upper weight instead of the p99 given by the proc univariate, used when p99 = 1 */
                  
           data_stats = 0 ,
		   seed = 1,
           debug = 0 /* save some of the intermediate files : alltrials, _tmp2_ */
		   
		   );


         %local m ;

         %if %bquote(&model_var) = %then %let model_var = &treatment ;

		 %if &check_for_cases = 0 and &check_controls = 1 %then %let final_analysis = 0 ;
		
         %if &data_stats = 1 %then %do;
              %let final_analysis = 0 ;
              %let check_for_cases = 0 ;
              %let check_controls = 0 ;
              %let checked_controls = 0 ;
              data for_stats ;
              run;
        %end;
 
        /* want to remove any included conditions in itt analysis. These two variables will be in the code when using itt analysis and censoring weights */
		%if &pp = 0 %then %do;
			%let condition0 = ;
			%let condition1 = ; 
		%end;	
		
        %if &run_user_limits_analysis = 1 and &run_p99_analysis = 1 %then %let run_p99_analysis = 0 ;

    	%let ncov = %numargs(&outcomeCov_var);
        data _datain_  %if &ncov > 0 %then (drop = _length_) ;;
        set &datain (keep = &id &time  &treatment &event &eligible_wts_0 &eligible_wts_1 
                            &cov_switchd &class_switchd
                           %if &eligible ^= _eligible_  %then &eligible ;
                           %if %bquote(&cense)^= %then &cense  &cov_censed &class_censed ;                         
                           %if %eval(&final_analysis) = 1   %then &outcomeCov_var  &where_var ;                        
                    )  ;

        %if &eligible = _eligible_ %then _eligible = 1 ;;
        if _n_ = 1 then do ;
          %do i = 1 %to &ncov ;
               %let word0 = %scan(&outcomeCov_var,%eval(&i),%str( ));
               _length_ = vlength(&word0);
               %global length&i ;
               call symput("length&i",trim(left(_length_)));
          %end;
       end; 
       run;

       

      
     /* assume that last observation has event = 1 or missing */ 



     proc sort data = _datain_ ;
	 by &id &time ;
	 run;



      data _for_outcome_;
      set _datain_ (where = (&event = 1));
      length _time_of_event 4 ;     
      _time_of_event = -1;
      if &event = 1 then _time_of_event = &time ;
      keep &id _time_of_event ;
      run;
	  
      data _datain_;
      merge _datain_    _for_outcome_;
      by &id;
      if _time_of_event = . then _time_of_event = -1 ;
      run;
         
        
      data _datain_  all_eligible_times (keep = &time);
      set _datain_ end = _end_;
      *retain  _regime_start   _cumA_ ; 
      by &id ;
     
      length switch   3 ;
                           
	  if &eligible = 1 then output all_eligible_times ;
      _Am1_ = lag(&treatment);
              
       if first.&id then do;                           
             _Am1_ = 0 ;
             switch = 0 ; /* set there to be no switch at "baseline" observation */                                 
       end;
       else do ;
            if _Am1_ ^= &treatment then switch = 1;
            else switch = 0 ;                            
       end;

                           
       *if _end_ then call symput("numids",trim(left(_newid_)));
	   output _datain_ ;
       run;   

     

	   proc sort data = all_eligible_times nodupkey ;
	   by &time ;
	   run;

	   proc sql noprint ;
	   select &time  into :all_eligible_times  separated by ' ' from all_eligible_times ;
       quit;
 

	   %put all m = &all_eligible_times ;
	   %let num_m = %numargs(&all_eligible_times );

	    %let first_eligible_m = %scan(&all_eligible_times , 1) ;
		%let last_eligible_m = %scan(&all_eligible_times,&num_m);


		%if &check_for_cases = 1 %then %do;
                %if &checked_for_cases = 0 %then %let m = %scan(&all_eligible_times,1);
                %else %if &checked_for_cases =1 %then %do;
				    proc sql ;
					select max(trial) as trial  into :m from &checked_data(keep = trial ) ;
					quit;

					%do i = 1 %to &num_m ;
					    %let testm = %scan(&all_eligible_times,&i);
						%if &m = &testm %then %let mcounter = &i ;
                    %end;
                    %put m = &m , place in list =  &mcounter  ;
				 %end;
	    %end ;


       %if &check_for_cases = 1 & &checked_controls = 1 %then  %let check_for_cases = 0 ;


		%if &check_for_cases = 1 and &final_analysis = 1 
              and (&m < &last_eligible_m or (&m = &last_eligible_m and &recheck_for_cases = 1))
			  %then %let final_analysis = 0 ;

	   /* when running after checking for cases _cases_ will contain any new cases that were found after cleaning the data 
	      this will be used for trial m+1 when recheck_for_cases = 0 */

	   proc sql noprint ;
	   create  table _cases_ as 
       select &id, &time ,&eligible,_time_of_event from _datain_ (  where = ( _time_of_event ne -1)) order by &time ;
	   quit;

	   

	   %let firstm = 1 ;
	   %let lastm = &num_m ;
       
       %let mtest = %scan(&all_eligible_times,1);
       %let mcases =0 ;
       %let i = 0 ;
       %do %while(&mcases = 0 AND &i < &lastm );
             %if &mcases =0 %then %do;
                %let i = %eval(&i + 1);
                %let mtest = %scan(&all_eligible_times,&i);
             %end;
       
            data _tmp0_ ;
            set _cases_ (where = (&time = &mtest and &eligible = 1));
            drop &eligible ;
            run; 

            %let mcases =%obscnt(_tmp0_) ;
          

       %end;
       
       %put first i with mcases > 0 is   &i ;
       %if &i = &lastm AND &mcases = 0 %then %do;
           %put NO CASES FOR ANY ELIGIBLE TIME ; 
           %goto exit ;
       %end;
       %else %do;
           %put RESETTING FIRSTM = &i FROM FIRST IN LIST OF ELIGIBLE TIMES DUE TO NO CASES IN FIRST ELIGIBLE TIME  ;
           %let firstm = &i ;
       %end;

	   %if &check_for_cases = 1 %then %do;
         %if &checked_for_cases = 0 %then %do;
		      %let firstm = 1 ;
			  %let lastm = 1;
		 %end;
		 %if &checked_for_cases = 1 %then %do;
              %let firstm = 1 ;
			  %if &recheck_for_cases = 1 %then %let lastm = 1 ;
			  %else %if &recheck_for_cases = 0 %then %let lastm = 2 ;
			  
			  %if &additional_controls = 0 %then %let firstm = 2 ; /* do not need any more controls. only need to add in next trial */
		 %end;

	   %end;

	   %if &checked_controls = 1 %then %let lastm = -1 ; /* want to skip the next loop since we have created all the trials */


    



    %do i = &firstm %to &lastm ;

        %if &check_for_cases =0 and &recheck_for_cases =0 %then %do;
            %let mm = %scan(&all_eligible_times,&i);
        %end;
        %if &check_for_cases = 1  %then %do;

            %if &i = 1 %then %let mm = &m ;
            %else %if &i = 2 %then %do;
                %let mcounter = %eval(&mcounter + 1);
                %let mm = %scan(&all_eligible_times,&mcounter ) ;
             %end;

        %end;		   


        %get_cases: %put determine cases for m = &mm ;

        %if &checked_for_cases = 0  or (&checked_for_cases = 1 and &mm > &m /* at a new eligible time */) %then %do;

            data _eligible_cases ;
            set _cases_ (where = (&time = &mm and &eligible = 1));
            drop &eligible ;
            run;

        %end ;
        %else %do ;
            data _eligible_cases ;
            set &checked_data (where = (trial = &mm and &event = 1)) ;
            run;
        


            /* time in this data is the followup time in the trial , qm is the selected time and corresponds to the the original
            data */
            data oldcontrols ;
            set &checked_data (keep = &id &time &event qm trial where = (trial = &mm and &event = 0)) ;
            keep &id qm ;
            run;

            proc sort data = oldcontrols ;
            by &id qm ;
            run;


        %end;

        %let mcases = %obscnt(_eligible_cases) ;
        /****/
        %if &check_for_cases = 1 AND &mcases = 0 AND &mm < &last_eligible_m %then %do;
            %let mcounter = %eval(&mcounter + 1);
            %let mm = %scan(&all_eligible_times,&mcounter ) ;
            %goto get_cases ;
        %end;
        /***/

        %if &mcases > 0 OR &data_stats = 1 %then %do;
            data _eligible_control_pool ;
            set _datain_ %if %bquote(&eligible_controls ) = %then (keep = &id &time &eligible &event where = (&time >= &mm  )) ;
			              %else (where = (time >= &mm)) ;;
            by &id ;
            retain _keep_;
            if first.&id then _keep_ = 0 ;
            if &time = &mm and &eligible = 1 then _keep_ = 1 ;
            if &event = 1 or &event = . then delete;
            if _keep_ = 0 then delete ;
            drop _keep_ &event &eligible ;
            run;
			
			%if %bquote(&eligible_controls )^= %then %do;
			
				data _eligible_control_pool ;
				set _eligible_control_pool ;
				if &eligible_controls  ;
				keep &id &time ;				
				run;
			
			
			%end;
			
			
            %let mcontrols = %obscnt(_eligible_control_pool);
            %if &data_stats = 1 %then %do;
                data for_stats ;
                set for_stats end = eof ;
                output ;
                if eof then do ;
                    m = &mm ; mcases = &mcases ; mcontrols = &mcontrols ; output ;
                end;
                run;

            %end;
            %if &checked_for_cases = 1 and &mm = &m %then %do; 
                /* want to remove the controls that we have already selected from the pool of
                    eligible ones */

/* this merge should add in qm for any id that was selected as a control. Will delete the observations where
time = qm since these observations were selected in a previous run of the macro */
                data _eligible_control_pool ;
                merge _eligible_control_pool (in = _a) oldcontrols(in = _b rename = (qm = &time)) ;
                by &id &time ;
                if _a = 1 and _b =1 then delete ;					
                run;

    
            %end;

            %let mcontrols =  %obscnt(_eligible_control_pool);

            %if &data_stats = 0 AND  &mcontrols >= %eval(&mcases * &ccratio) %then %do;
                /* using srs  will select at most one time, numberhits = 1 for those selected, can remove numberhits from robust macro */
                proc surveyselect data = _eligible_control_pool noprint
                method = srs 
                seed = %eval(&seed * &mm + 1)
                %if &checked_for_cases = 0 or (&checked_for_cases and &mm > &m /* at a new eligible time */)  %then n = %eval(&mcases * &ccratio) ;
                %if &checked_for_cases = 1 and &mm = &m %then n = &additional_controls ;
                out = selected_controls (rename = (&time = qm)  );
                run;

                data trialm ;
                %if &checked_for_cases =  0 or (&checked_for_cases = 1 and (&mm > &m)) %then set _eligible_cases selected_controls ;;
                %if &checked_for_cases = 1 and (&mm = &m) %then 
                set 
                    &checked_data (keep = &id &time qm &event trial  where = (trial = &mm )) /* previously selected cases and controls */ 
                    selected_controls /* additional controls */ ;;
                trial = &mm ;	  			
                %if &checked_for_cases =  0 or (&checked_for_cases = 1 and (&mm > &m)) %then %do; 
                    if _time_of_event > . then &event = 1;
                    else &event = 0 ;
                    if &event = 1 then &time = _time_of_event - trial ;
                    else if &event = 0 then &time = qm - trial ;

                %end;
                %if &checked_for_cases = 1 and (&mm = &m) %then  %do;
                    if &event = . then &event = 0 ; /* for added controls with no event variable */
                    &time = qm - trial ;
                %end;
                if &event = 1 then qm = &time + trial ;
                %if &checked_for_cases =  0 or (&checked_for_cases = 1 and (&mm > &m)) %then  drop _time_of_event ;; 
                run;

                proc sort data = trialm ;
                by &id ;
                run;


                %if &check_for_cases = 0 %then %do;
                    %if &i = &firstm %then %do;
                        data alltrials ;
                        set trialm ;
                        run;
                    %end;
                    %else %do;
                        proc append base = alltrials data = trialm ;
                        run;
                    %end;
                %end;
                %else %if &check_for_cases = 1 %then %do;
                        /* this should work for any pass through macro, including the first */
                    %if &mm = &first_eligible_m %then %do;
                        data alltrials ;
                        set trialm ;
                        run;
                    %end;
                    %if &mm > &first_eligible_m %then %do;
                        %if &mm = &m %then %do; /* first iteration of loop: could be first try for trial m or for additional recheck */
                            data alltrials ;
                            set &checked_data (where = (trial < &m)) /* data from previous trials */
                                trialm  /* updated data for trial &m */   ; 
                            run;
                        %end;
                        %else %if &mm > &m %then %do ; /* for a possible second iteration of loop for new trial */
                            %if &additional_controls > 0 %then %do;
                                proc append base = alltrials data = trialm ;
                                run;
                            %end;
                            %else %if &additional_controls = 0 %then%do;
                                data alltrials ;
                                set &checked_data /* use whole data set since we did not need to add in any new controls */ 
                                    trialm ;
                                run;
                            %end;

                        %end; 
                    %end;

                %end;
            %end;
        %end;  


    %end;

    %if &debug = 1 %then %do;
           data initial_alltrials ;
           set alltrials ;
           run;
    %end;
 

%if &final_analysis = 1 %then %do; 


        %if  &checked_controls = 1 %then %do;
            %let check_for_cases = 0 ;
            proc sort data = &checked_data out = _tmp_ ;
			by &id trial ;
			run;

		
		%end;
		%else %do;

            proc sort data = alltrials out = _tmp_ ;
			by &id trial ;
			run;
 

		 

		%end;

        %if &pp = 0 AND %bquote(&cense) =  %then %do;
            /* for itt analysis, add in the covariates for each case and control */
	       data alltrials ;
  		    merge _tmp_ (in = __a) _datain_ ( drop = &event switch _Am1_ rename = (&time = trial ));
      	    by &id trial ;
  		    if __a ;
  			run;

        %end;
        %else %if &pp = 1 OR %bquote(&cense)^=  %then %do;

            %let min_eligible = &first_eligible_m ;
            * this data creation only saves cases if they are selected as controls at least once ;
            proc sql ;
            create table for_controls as 
            select &id, min(trial) as min_selected0 ,max(qm) as max_selected from _tmp_(where =(&event =0))  group by &id ;
            create table for_cases as 
              select &id, min(trial) as min_selected1 ,max(qm) as max_selected1 from _tmp_(where =(&event =1))  group by &id ;
            quit;



            data max_selected  ;
            merge for_controls for_cases ;
            by &id ;
            changed = 0 ;
           
            *if min_selected = . then min_selected = min_selected1 ;
            min_selected = min(min_selected0,min_selected1) ;

            /* No observations selected as control, so all observations  are  for cases. 
               In data for models should have the outcome missing  */
            if max_selected = . then do ; 
                ***max_selected = max_selected1 ;                
                   max_selected = -1 ; /* assumed minimum time in data is at least 0 */                
            end;

            keep &id min_selected max_selected ;
            run;
 
            %if &pp = 1 %then %do; /* will use treatment models */

                data _tmp2_ ;
                merge _datain_ max_selected(in = _a) ;
                by &id ;
                if _a = 1 ;						
                _current_treatment = &treatment ;


             

                %if %bquote(&cense) = %then %do;
		              /* for treatment weights and no censoring, keep _tmp2_ the way it was */
				        * keep observations between min(m) and max_selected for controls and between min(m) and event time for cases ;
				     if (_time_of_event = -1 and  min_selected <= &time <= max_selected ) or 
					   (_time_of_event ne -1 and min_selected  <= &time <= _time_of_event) ;
				    *** if _time_of_event ne -1 then &treatment = . ;
				    * for cases we want to set the treatment to be missing for times after the last selected time point and the time of the event ;
				    if _time_of_event ne -1  and max_selected < &time <= _time_of_event then  &treatment = . ;
                %end;
			    %else %do ; /* for both treatment and censoring models */
				
			         if min_selected <= &time ; /* include all observations after minimum selected m  */
					
				     if _time_of_event = -1 and  max_selected < &time then &treatment = . ; /* do not want to include observations after max qm in treatment models */				
        		     if _time_of_event ne -1  and max_selected< &time <= _time_of_event then do;
					     &treatment = .; /* remove from treatment models */
						 &cense = . ;   /* remove from censoring models */						 
				    end;				             
			   %end;							
               run;
           
            
            %if  &pool_models = 0 %then %do; 
                title "Treatment model given lagged value is 0 : denominator";
                proc logistic data = _tmp2_ (where = (_Am1_ = 0)) ;
			    %if %bquote(&class_switchd)^= %then class &class_switchd   ;;
				ods select ModelInfo ResponseProfile ParameterEstimates ;
                model &treatment (descending) = &cov_switchd ;
                output out=_weights0d_   pred=pA_d ;
                run;
                title "Treatment model given lagged value is 1 : denominator";
                proc logistic data = _tmp2_ (where = (_Am1_ = 1)) ;
				ods select ModelInfo ResponseProfile ParameterEstimates ;
			    %if %bquote(&class_switchd)^= %then class &class_switchd   ;;
                model &treatment (descending) = &cov_switchd ;
                output out=_weights1d_   pred=pA_d ;
                run;
                 title ;
            
				%if &use_stabilized_weights = 1  %then %do;
				
					title "Treatment model given lagged value is 0 : numerator";
					proc logistic data = _tmp2_ (where = (_Am1_ = 0)) ;
					ods select ModelInfo ResponseProfile ParameterEstimates ;
					%if %bquote(&class_switchn)^= %then class &class_switchn   ;;
					model &treatment (descending) = &cov_switchn ;
					output out=_weights0n_   pred=pA_n ;
					run;
					title "Treatment model given lagged value is 1 : numerator";
					proc logistic data = _tmp2_ (where = (_Am1_ = 1)) ;
					ods select ModelInfo ResponseProfile ParameterEstimates ;
					%if %bquote(&class_switchn)^= %then class &class_switchn   ;;
					model &treatment (descending) = &cov_switchn ;
					output out=_weights1n_   pred=pA_n ;
					run;
					title ;
				
					data _weightsn_;
					set _weights0n_ _weights1n_;
					by &id &time;
				
				%end;
			
			
                data _weights_ ;
                set _weights0d_ _weights1d_ ;
                by &id &time ;
                run;
				
				%if &use_stabilized_weights = 1 %then %do;
					data _weights_ ;
					merge _weights_ _weightsn_;
					by &id &time ;
					run;
				
				%end;
            %end ;
            %else %if &pool_models = 1 %then %do;
                title "Treatment model pooled over previous values being 0 and 1 : denominator ";
                proc logistic data = _tmp2_   ;
				ods select ModelInfo ResponseProfile ParameterEstimates ;
			    %if %bquote(&class_switchd)^= %then class &class_switchd   ;;
                model &treatment (descending) = &cov_switchd ;
                output out=_weights_   pred=pA_d ;
                run;
                 title ;
				 
				 %if &use_stabilized_weights = 1 %then %do;
					title "Treatment model pooled over previous values being 0 and 1 : numerator ";
					proc logistic data = _tmp2_   ;
					ods select ModelInfo ResponseProfile ParameterEstimates ;
					%if %bquote(&class_switchn)^= %then class &class_switchn   ;;
					model &treatment (descending) = &cov_switchn ;
					output out=_weightsn_   pred=pA_n ;
					run;
					title ;
					
					
					data _weights_;
					merge _weights_ _weightsn_ ;
					by &id &time ;
					run;
				 
				 %end;

            %end;

            %if &debug = 1 %then %do;
                data treatment_model_data ;
                set _tmp2_;
                run;
            %end;

			%end;
			



			%if %bquote(&cense) ^= %then %do ;
                %if &pp = 0 %then %do;
                    data _tmp2_ ;
                    merge _datain_ max_selected(in = _a) ;
                    by &id ;
                    if _a = 1 ;						
                    _current_treatment = &treatment ;
				
			         if min_selected <= &time ; /* include all observations after minimum selected m  */
					
				    * make no changes to the true controls in the data, use complete follow-up from first selected m ;
                     * for controls from cases, only use observations untill the maximum selected qm for censoring model ;
        		     if _time_of_event ne -1  and max_selected< &time <= _time_of_event then do;
					    
						 &cense = . ;   /* remove from censoring models */						 
				    end;				             
			  							
                    run;
           


                %end;
                %else %do;
                    *need to reset the treatment variable for the controls from the cases, can do this for all observations ;

                    data _tmp2_;
                    set _tmp2_ ;
                    &treatment = _current_treatment ;
                    run;
                 %end;
				proc logistic data =  _tmp2_      ;  
				ods select ModelInfo ResponseProfile ParameterEstimates ;                                 
				title3 "Model for P(&cense = 0 |  X ) for denominator ";
				class  _Am1_ ( ref = first )  &class_censed / param=ref;
				model &cense (ref = last)  = _Am1_ &cov_censed          ;                          
				output out =  _cense_ %if &pp = 1 %then (keep = &id &time pC_d &cense )  ;      PRED = pC_d ;
				run; 



				%if &use_stabilized_weights  = 1 %then %do;
					proc logistic data =  _tmp2_    ; 
					ods select ModelInfo ResponseProfile ParameterEstimates ;
					title3 "Model for P(&cense = 0 |  X ) for numerator ";
					class _Am1_ ( ref = first)  &class_censen  /param=ref;
					model &cense (ref = last) = _Am1_ &cov_censen        ;                                
					output out =  _cense_n0 (keep = &id &time pC_n &cense )       PRED = pC_n ;
					run; 
				
					data _cense_ ;
					merge _cense_ _cense_n0 ;
					by &id &time ;
					keep &id &time pC_d pC_n %if &pp = 0 %then _current_treatment; ;
					if &cense = 1 then do ;
						pC_d = 1 - pc_d ;
						pc_n = 1 - pC_n ;
					end;	
					run;
				 
				
				%end ;
		   
				%if &pp = 1 %then %do;
					data _weights_ ;
					merge _weights_ _cense_ ;		  
					by &id &time ;
					run;
				%end;
				%else %do;
					data _weights_;
					set _cense_;
					run;
				%end;

                %if &debug = 1 %then %do;
                    data censoring_model_data ;
                    set _tmp2_;
                    run;
                %end; 
			%end;
			
			 %let first_eligible_m = %sysfunc(left(&first_eligible_m));
			
			%if %bquote(&cense)= %then %do;
                title 'maximum selected qm for all cases/controls';
                proc sql ;
                select max(qm) as max_selected into :max_selected_time from _tmp_ ;             
                quit;
                title ;
           
                %let max_selected_time = %sysfunc(left(&max_selected_time));
            %end;
            %else %do;
                 title 'for censoring/treatment weights maximum selected time is maximum time '; 
                 proc sql ;
                 select max(&time) as max_selected into :max_selected_time from _tmp2_ ;             
                 quit;
                 title ;
                %let max_selected_time = %sysfunc(left(&max_selected_time));
             
            %end;
            

            /* each row of _weights_ contains the history of _wt and _at for the previous times.
               we can use this,as in the initiators macro, to calculate the weight for weight at qm for trial m 
               and to determine if there was a change in treatment between time = m and time = qm.  */
            data _weights_  ;
            set _weights_ ;
            by &id ;
            retain _wt&first_eligible_m -_wt&max_selected_time  _at&first_eligible_m -_at&max_selected_time ; 
            array _wt{&first_eligible_m:&max_selected_time} _wt&first_eligible_m -_wt&max_selected_time ;
            array _at{&first_eligible_m:&max_selected_time} _at&first_eligible_m -_at&max_selected_time ;
            %if %bquote(&condition0)^= %then %do;
                  array _condition0_{&first_eligible_m:&max_selected_time} _condition0_&first_eligible_m -_condition0_&max_selected_time ;
                  retain _condition0_&first_eligible_m -_condition0_&max_selected_time ;
                  if first.&id then do ;
                    do i = &first_eligible_m to &max_selected_time ;
                       _condition0_[i] = 0 ;                      
                    end;
                  end;
                  _condition0_[&time] = (&condition0 ) ;
                  
            %end;
                    
            %if %bquote(&condition1)^= %then %do;
                  array _condition1_{&first_eligible_m:&max_selected_time} _condition1_&first_eligible_m -_condition1_&max_selected_time ;
                  retain _condition1_&first_eligible_m -_condition1_&max_selected_time ;
                  if first.&id then do ;
                    do i = &first_eligible_m to &max_selected_time ;
                       _condition1_[i] = 0 ;                      
                    end;
                  end;
                  _condition1_[&time] = (&condition1 ) ;
            %end; 
            if first.&id then do ;
                do i = &first_eligible_m to &max_selected_time ;
                     _wt[i] = 1.0 ;
                     _at[i] = . ;
                end;
            end;
            _at[&time] = _current_treatment ;
            * for times between max sected time and time of event, &treatment is missing. will need to use _current_treatment instead ;
			numerator = 1.0 ;
			denominator = 1.0 ;
			 %if &pp = 1 %then %do;
				if _current_treatment  = 0 then pA_d = 1.0 - pA_d ;
				denominator = pA_d ;
			%end;	
			
			%if %bquote(&cense)^= %then %do;
			     denominator = denominator * pC_d ;
			%end;	 
           
			%if &use_stabilized_weights = 1 %then %do;
				%if &pp = 1 %then %do;
					if _current_treatment = 0 then pA_n = 1.0 - pA_n ;
					numerator = pA_n ;
				%end;	
				%if %bquote(&cense)^= %then %do;
				    numerator = numerator * pC_n ;
				%end;				    		
			%end;	
			 
            if &time = &first_eligible_m then _wt[&time] = numerator/denominator ;
            else if &time > &first_eligible_m then _wt[&time] = _wt[&time - 1] * (numerator/denominator) ;
            *drop i &condition_vars &treatment _current_treatment ;
            keep &id &time  _wt&first_eligible_m -_wt&max_selected_time 
                 _at&first_eligible_m -_at&max_selected_time  
                 %if %bquote(&condition0)^= %then   _condition0_&first_eligible_m -_condition0_&max_selected_time ;
                 %if %bquote(&condition1)^= %then   _condition1_&first_eligible_m -_condition1_&max_selected_time ;
                 ;
            run;

/* _tmp3_ is a reduced version of _weights_ to contain only the selected time points for all the trials */            
proc sql ;
create table _tmp3_ as select 
 a.*  , b.trial, b.&time , b.&event 
   from _weights_ (rename = (&time = qm) )  as a right join _tmp_ as b 
   on( a.&id = b.&id and  a.qm = b.qm )
    order by &id , trial  ;
   quit;

 /* use _at to determine if there is a change in treatment from time = m to time = qm for each m, 
    this will be contained in the anyswitch variable and will be used for censoring subjects in trial m */
	
	data _tmp4_ ;
	set _tmp3_ ;
	array _wt{&first_eligible_m:&max_selected_time} _wt&first_eligible_m -_wt&max_selected_time ;
	_wt_ = _wt[qm]/_wt[trial] ;
	drop _wt&first_eligible_m -_wt&max_selected_time ;
	
	array _at{&first_eligible_m:&max_selected_time} _at&first_eligible_m -_at&max_selected_time ;
	&treatment = _at[trial] ;
	
	%if &pp = 1 %then %do;
		
		anyswitch = 0 ;
		time_of_switch = -1 ;
		do i = trial + 1 to qm ;
			if _at[i] ne _at[trial] and anyswitch = 0 then do;
				anyswitch = 1 ;
				time_of_switch = i ;
			end;
		end; 
		
	%end;
	drop _at&first_eligible_m -_at&max_selected_time  ;
	
	%if %bquote(&condition0)^= %then %do;
          
          array _condition0_{&first_eligible_m:&max_selected_time} _condition0_&first_eligible_m -_condition0_&max_selected_time ;
          _anycond0 = 0;
          _time_of_cond0 = -1 ;
          if &treatment = 0 then do ; * treatment is value at start of trial , time = m = trial ;
            do i = trial to qm ;
                  if _condition0_[i]= 1 and _anycond0 = 0  then do ;
                        _anycond0 = 1 ;
                        _time_of_cond0 = i ;
                   end;
            end;
          end;
          if _anycond0 = 1 and (time_of_switch >= _time_of_cond0) then do;
               anyswitch = 0 ;
                _wt_ = _wt[_time_of_cond0]/_wt[trial]; *redefine wt to be wt at time of condition 0 satisfied ;
         
          end;
          drop _condition0_&first_eligible_m -_condition0_&max_selected_time _anycond0 _time_of_cond0 ;         
	%end;
	%if %bquote(&condition1)^= %then %do;
          
          array _condition1_{&first_eligible_m:&max_selected_time} _condition1_&first_eligible_m -_condition1_&max_selected_time ;
          _anycond1 = 0;
          _time_of_cond1 = -1 ;
          if &treatment = 1 then do ; * treatment is value at start of trial , time = m = trial ;
            do i = trial to qm ;
                  if _condition1_[i]= 1 and _anycond1 = 0  then do ;
                        _anycond1 = 1 ;
                        _time_of_cond1 = i ;
                   end;
            end;
          end;
          if _anycond1 = 1 and (time_of_switch >= _time_of_cond1) then do;
               anyswitch = 0 ;
               _wt_ = _wt[_time_of_cond1]/_wt[trial]; *redefine wt to be wt at time of condition 1 satisfied ;
         
          end;
          drop _condition1_&first_eligible_m -_condition1_&max_selected_time _anycond1 _time_of_cond1 ;         
	%end;
	/* drop i  ; */
	run;

	%if &pp = 1 %then %do;
		title "number of observations censored from the pp data set due to changing treatment";
		proc freq data= _tmp4_ ;
		table anyswitch ;
		run;
		title ;

        %if &debug = 1 %then %do;
             data pre_switch_censoring ;
             set _tmp4_;
             run;

        %end;

		data _tmp4_;
		set _tmp4_ (keep = &id _wt_ anyswitch  &time trial qm &event &treatment ) ;
		if anyswitch = 1 then delete ; /* censor those who change treatment in trial = m */
		run; 
		
		
	%end;
		
    /* when merged treatment was set from arrays in _tmp4_ */
	data alltrials ;
	merge _tmp4_ (in = __a) _datain_ ( drop =   &treatment  &event switch _Am1_ rename = (&time = trial    ));
	by &id trial ;
	if __a ;
	run;

	proc means data = alltrials min p1 p10 p25 p50 p75 p90 p99 max mean std ;
	var _wt_ ;
	output out= p (keep = p1 p99) p1=p1 p99=p99 ;
	run;

     %if &run_p99_analysis = 1 %then %do; 
        data _null_ ;
        set p ;
        call symput('p01',compress(p1));
        call symput('p99',compress(p99));
        run;

        data alltrials ;
        set alltrials(rename = (_wt_ = _wtorig_));
        _wt_ = _wtorig_;
        if _wt_ >= &p99 then _wt_ = &p99 ;
        else if _wt_ <= &p01 then _wt_ = &p01;
        run ;


        proc means data = alltrials min p1 p10 p25 p50 p75 p90 p99 max mean std ;
        var _wt_ ;
  
        run;
     %end;

     %if &run_user_limits_analysis = 1 %then %do; 
         

        data alltrials ;
        set alltrials(rename = (_wt_ = _wtorig_));
        _wt_ = _wtorig_;
        if _wt_ >= &upper_weight then _wt_ = &upper_weight ;
        else if _wt_ <= &lower_weight then _wt_ = &lower_weight;
        run ;


        proc means data = alltrials min p1 p10 p25 p50 p75 p90 p99 max mean std ;
        var _wt_ ;
  
        run;
     %end;


%end ;

	%let drop_list =  &event   &model_var  trial &time   &outcomeCov_var  qm  ;
    %if &pp = 1 OR %bquote(&cense)^= %then %let drop_list = &drop_list _wt_ ;

    %let keep_list = &id &time trial &treatment &event &outcomecov qm ;
    %if &pp = 1 OR %bquote(&cense)^=  %then %let keep_list = &keep_list _wt_ ;

	proc logistic data = alltrials (keep = &keep_list ) descending   outest = _est_block covout ;
	 %if %bquote(&outcomeClass)^= %then class &outcomeClass   ;;	                   	              
	 model  &event =  &model_var                      
	                 %if %eval(&include_followup_time_outcome) = 1 %then &time ; 
	                 %if &include_expansion_time_outcome = 1 %then trial ;
	                 &outcomeCov ; 

     %if &pp = 1 OR %bquote(&cense)^= %then weight _wt_ ;;
      %if &calculate_var = 1 %then %do;
            output out =_for_robust(   drop = &drop_list )   dfbeta = _ALL_  ;
      %end;
     run; 


	              data _est_tmp;
	              set _est_block ;
	              if _n_ = 1 ;                            
	              drop _LINK_ _TYPE_ _STATUS_ _NAME_ _LNLIKE_ _ESTTYPE_;
	              run;

	              %if &calculate_var = 1 %then %do;                        
	 
	                   data _for_dfbetas;
	                   set _for_robust(obs = 1);
	                   run;

	                   proc transpose data = _for_dfbetas out = _for_dfbetast;
	                   run;

	                   data _for_dfbetast;
	                   set _for_dfbetast (keep = _NAME_);
	                   if lowcase(_NAME_) = lowcase("&id") then delete;
					   if lowcase(_NAME_) = "numberhits" then delete ;
	                   run;

	                   proc sql noprint;
	                   select _NAME_
	                   into :for_dbeta separated by ' '
	                   from _for_dfbetast;
	                   quit;

	                    %robust(outest=_est_block, out=_for_robust, id=&id, df=&for_dbeta);
	               
	                 /* output of robust macro is a data set named se_block and covar_block */

	                    proc datasets library = work nolist;
	                    delete   _for_robust _est_block _reduce_ _out1_  ;
	                    quit;
	                   
						    
	   
	                           data _tmp;
	                           set _est_tmp;
	                           run;

	                           proc transpose data = _tmp out= _tmpt /*(drop = COL1) */;
	                           run;

	                           data _tmpt ;
	                           set _tmpt;
	                           if _LABEL_ = '' then _LABEL_ = _NAME_;
	                           if _n_ = 1 then _LABEL_ = 'intercept';
							   _var_ = _N_ ;
	                           run;
	       
	                           proc transpose data = se_block out =se_t (drop = _NAME_) ;
	                           run;

							   data se_t ;
							   set se_t ;
							   _var_ = _N_;
							   run;

	                          data output    ;
	                          merge  _tmpt ( rename = (col1 = estimate _NAME_ = parameter)) se_t (rename = (col1 = std));
	                         
	                          by _var_ ;
	                          label p_value = 'Pr > |Z|';
	                          format p_value PVALUE8.4 ;  
	                          if std > 0 then do ;   
	                              lb = estimate - 1.96 * std;
	                              ub = estimate + 1.96 *std ;
	                              z = estimate / std;
	                              p_value = 2*(1-probnorm(abs(z))); 
	                          end;
	                          else std = . ;
							  drop _var_ ;
	                          run;
	                       %if &pp=0 AND %bquote(&cense)= %then  title "Analysis of ITT with no weights and using robust variance";;
                           %if &pp=0 AND %bquote(&cense)^= %then  title "Analysis of ITT with censoring weights and using robust variance";;
                           %if &pp=1 AND %bquote(&cense)= %then title "Analysis of ITT with treatment weights and using robust variance";;
                           %if &pp=1 AND %bquote(&cense)^= %then title "Analysis of ITT with treatment and censoring weights and using robust variance";;
	                       proc print data = output noobs;
	                       var parameter estimate std lb ub z p_value ;
	                       run;
						   title ;

	        
	               %end;  
	       
%end;
%if &data_stats = 1 %then %do;
    data for_stats ;
    set for_stats ;
    if m ne . ;
    label m='Eligible time (m)'
          mcases = 'Number of cases'
          mcontrols = 'Size of control pool';
    run;

    proc print data = for_stats noobs label;
    var m mcases mcontrols ;
    run ;

%end;

%exit: ;

%mend ;
%macro robust(outest=a, out=_last_, id=id, df= );

     /*********************************************************************
     **********************************************************************
     MACRO ROBUST produces robust estimates of the covariance matrix for
     the coefficients from PROC LOGISTIC or PROC PHREG when the data are
     clustered. One common application is to longitidunal data with each
     individual treated as a cluster.  IML is required.

     The macro uses the method of Halbert White (1982) "Maximum likelihood
     estimation of misspecified models," Econometrica 50: 1-25.

     Author:  Paul D. Allison, University of Pennsylvania
     allison@ssc.upenn.edu

     Adapted from an IML program written by Terry Therneau, Mayo Clinic.

     For PROC LOGISTIC, you must specify the OUTEST=name1 and
     COVOUT options on the PROC statement. You must also use the OUTPUT
     statement with OUT=name2 and DFBETAS=namelist.  The namelist should
     have one name for each term in the model, including the intercept.
     There must also be a variable in the data set containing a unique
     value (either character or numeric) for each cluster. The macro is
     invoked after running the LOGISTIC procedure.

     For PROC PHREG, you must specify OUTEST=name1 on the PROC
     statement. (COVOUT is unncessary).  You must also use the OUTPUT
     statement with OUT=name2 and DFBETA=namelist.  The namelist should
     have one name for each variable in the model. There must also be a
     variable in the data set containing a unique value (either character
     or numeric) for each cluster. This variable must be added to the
     OUTPUT data set by using the ID statement.  The macro is
     invoked after running the PHREG procedure.

     The macro has the following parameters:

     OUTEST   Name of data set used in the OUTEST= option.
     OUT      Name of data set used in the OUT= option.
     ID       Name of variable containing a unique value for each cluster
     DF       List of names used in the DFBETAS or DFBETA option.

     Examples of usage:

     proc logistic outest=a covout;
     model y = x z w;
     output out=b dfbetas=dint dx dz dw;
     run;
     %robust(outest=a, out=b, id=subjid, df=dint dx dz dw)

     proc phreg outest=a;
     model y*d(0) = x z w;
     id subjid;
     output out=b dfbeta= dx dz dw;
     run;
     %robust(outest=a, out=b, id=subjid, df=dx dz dw)

     BE CAREFUL: it's DFBETAS in LOGISTIC but DFBETA in PHREG.
     (The former is standardized, the latter is not).
     Also PHREG does NOT have an intercept.
     **************************************************************
     **************************************************************/


     proc means data=&out noprint;
     class &id;
     var &df;
     output out=_out1_(keep=&df) sum=&df;
     types &id ;
	 
     run;


     data _out1_;
     set _out1_ ;
     array d (*) &df;
     if d(1)=. then delete;
     do _i_ = 1 to dim(d) ;
         if d(_i_) = . then d(_i_) = 0 ;
     end;
     drop _i_ ;
     run;

     data _reduce_;
     set &outest;
     array abc(*) _character_;
     length name $8;
     call vname(abc(1),name);
     if name ne '_LINK_' and _type_ eq 'COV' then delete;
     drop  _lnlike_;
     run;


     proc iml;
     use _reduce_ where (_type_='COV');
     read all into cov;
     use _reduce_ where (_type_='PARMS');
     read all into b[colname=vname];
     if ncol(cov)=0 then se=1;
     else se=sqrt(diag(cov));

     nn = ncol(se);
     do i = 1 to nn ;
        if se[i,i] = . then se[i,i] = 0 ;
     end;

     use _out1_;
     read all into x;
     x=x*se;
     v=x`*x;
     se=sqrt(vecdiag(v));
     reset noname fuzz=.000001;

     se = t(se);
     create se_block from se ;
     append from se ;

     create covar_block from v;
     append from v ;

     quit;
     run;
   

%mend robust;
%macro obscnt(dsn);
%local nobs dsnid;
%let nobs=.;
 
%* Open the data set of interest;
%let dsnid = %sysfunc(open(&dsn));
 
%* If the open was successful get the;
%* number of observations and CLOSE &dsn;
%if &dsnid %then %do;
     %let nobs=%sysfunc(attrn(&dsnid,nlobs));
     %let rc  =%sysfunc(close(&dsnid));
%end;
%else %do;
     %put Unable to open &dsn - %sysfunc(sysmsg());
%end;
 
%* Return the number of observations;
&nobs
%mend obscnt;

%macro numargs(arg);

     %let n = 1;
     %if &arg^= %then %do;
          %do %until (%scan(&arg,%eval(&n),%str( ))=%str());
               %let word = %scan(&arg,&n);
               %let n = %eval(&n+1);
          %end;
     %end;
     %eval(&n-1) /* there is no ; here since it will be used as %let a = %numargs(&b) ;
     and the ; is included at the end of this line  */                            
                            
%mend numargs;
