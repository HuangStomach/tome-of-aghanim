package com.example.multidatasource.config;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.autoconfigure.orm.jpa.JpaProperties;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.jdbc.DataSourceBuilder;
import org.springframework.boot.orm.jpa.EntityManagerFactoryBuilder;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Primary;
import org.springframework.data.jpa.repository.config.EnableJpaRepositories;
import org.springframework.orm.jpa.JpaTransactionManager;
import org.springframework.orm.jpa.LocalContainerEntityManagerFactoryBean;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;

import javax.persistence.EntityManagerFactory;
import javax.sql.DataSource;
import java.util.HashMap;
import java.util.Map;

@Configuration
@EnableJpaRepositories(
    basePackages = "com.example.multidatasource.repository.a",
    entityManagerFactoryRef = "entityManagerFactoryA",
    transactionManagerRef = "transactionManagerA"
)
public class DataSourceConfigA {
    @Primary
    @Bean(name = "dataSourceA")
    @ConfigurationProperties(prefix = "spring.datasource.a")
    public DataSource dataSourceA() {
        return DataSourceBuilder.create().build();
    }

    @Bean(name = "jpaPropertiesA")
    @ConfigurationProperties(prefix = "spring.jpa.a")
    public JpaProperties jpaPropertiesA() {
        return new JpaProperties();
    }

    @Primary
    @Bean(name = "entityManagerFactoryA")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryA(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceA") DataSource dataSource) {
        Map<String, String> props = jpaPropertiesA().getProperties();

        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.entity.a")
                .persistenceUnit("a")
                .properties(props)
                .build();
    }

    @Primary
    @Bean(name = "transactionManagerA")
    public PlatformTransactionManager transactionManagerA(
            @Qualifier("entityManagerFactoryA") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }
}
